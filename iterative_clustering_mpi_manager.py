from mpi4py import MPI
import queue
import pickle
import time
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc

from transcriptomic_clustering.iterative_clustering import onestep_clust, OnestepKwargs

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MANAGER_RANK = 0

def setup_transcriptomic_clustering(): 
    means_vars_kwargs = {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, 0.6931472 to match to bigcat
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    }
    highly_variable_kwargs = {
        'max_genes': 4000 # originally 3000, 4000 to match to bigcat
    }
    pca_kwargs = {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    }
    ## Leave empty if you don't want to use known_modes
    filter_known_modes_kwargs = {
        # 'known_modes': known_modes_df, # a pd dataframe. index is obs (cell) names, columns are known modes. Originally commented out 
        'similarity_threshold': 0.7
    }
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    latent_kwargs = {
        # 'latent_component': "X_pca"
        'latent_component': "scVI"
    }

    cluster_louvain_kwargs = {
        'k': 15, # number of nn, originally 150, change to 15
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus # originally 8`
        'resolution': 1.0 # resolution of louvain for taynaud method
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, 10 to match to bigcat
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 1, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, 100 to match to bigcat
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, 0.6931472 to match to bigcat
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, 4 to match to bigcat
        'de_method': 'ebayes'
    }

    onestep_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )
    return onestep_kwargs

def manager_job_queue(adata_path, scvi_path, out_path): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()

    adata = sc.read(adata_path)
    print(f"Finished reading in anndata: {adata}")

    if np.max(adata.X) > 100:
        print(f"raw count data provided")
        print(f"Normlazing total counts to 1e6...")
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        print(f"Finished normalization. max:{np.max(adata.X)}")
    else:
        print(f"normalized data provided")

    latent = pd.read_csv(scvi_path, index_col=0)
    latent = latent.loc[adata.obs_names]
    print(f"Finished reading in scVI latent space: {latent.shape}")
    
    adata.obsm['scVI'] = np.asarray(latent)

    kwargs = setup_transcriptomic_clustering()
    min_samples = 4
    random_seed = 2024
    min_samples = 4
    tmp_dir =  os.path.join(out_path,'tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # let the manager do the first clustering
    clusters, markers = onestep_clust(adata, kwargs, random_seed)
    sizes = [len(cluster) for cluster in clusters]
    print(f"Manager finished clustering {adata.shape[0]} cells into {len(clusters)} clusters of sizes {sizes}")

    results = []
    
    # Initialize job queue, and put all subclusters into the queue
    adata_tmp_idx = 1
    job_queue = queue.Queue()
    if(len(clusters) == 1):
        results.append(clusters[0]) 
    else:
        for i in range(len(clusters)):
            if len(clusters[i]) < min_samples:
                results.append(clusters[i])
            else:
                idx = clusters[i]
                new_adata = adata[idx]
                new_adata_path = os.path.join(tmp_dir, str(adata_tmp_idx)+'.h5ad')
                new_adata.write(new_adata_path)

                new_task = (onestep_clust, (new_adata_path, kwargs, random_seed, idx))
                job_queue.put(new_task)

                adata_tmp_idx += 1

    # send the clusters from the first one-step clustering for further clustering 
    # there might be clustering tasks left in the queue
    active_workers = 0

    for worker_rank in range(1, size):
        if not job_queue.empty():
            task = job_queue.get()
            comm.send(task, dest=worker_rank, tag=1)
            active_workers += 1
            # print(f"Master sent a clustering task to Worker {worker_rank}")
        else:
            comm.send(None, dest=worker_rank, tag=0) # terminate nodes that have no tasks
            # break
    print(f"Initiate {active_workers} activate workers")

    # manage completed tasks from workers and assign new tasks
    while not job_queue.empty() or active_workers > 0:
        status = MPI.Status() # it can be defined outside of the while loop, does not matter
        clusters, new_markers, n_cells = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status) # receive any signals from any slaves
        worker_rank = status.Get_source() # determine which worker just sent a message to the master
        # tag = status.Get_tag() # tag==1: actual results, tag==0: acknowledgement of termination signal

        markers = markers|new_markers

        sizes = [len(cluster) for cluster in clusters]
        print(f"Worker {worker_rank} finished clustering {n_cells} cells into {len(clusters)} clusters of sizes {sizes}")

        # save the finished clustering results and add new tasks of cells to be clustered
        if(len(clusters) == 1):
            results.append(clusters[0])
        else:
            for i in range(len(clusters)):
                if len(clusters[i]) < min_samples:
                    results.append(clusters[i])
                else:
                    idx = clusters[i]
                    new_adata = adata[idx]
                    new_adata_path = os.path.join(tmp_dir, str(adata_tmp_idx)+'.h5ad')
                    new_adata.write(new_adata_path)

                    new_task = (onestep_clust, (new_adata_path, kwargs, random_seed, idx))
                    job_queue.put(new_task)

                    adata_tmp_idx += 1

        if not job_queue.empty():
            new_task = job_queue.get()
            comm.send(new_task, dest=worker_rank, tag=1)
        else:
            active_workers -= 1
            comm.send(None, dest=worker_rank, tag=0)
            print(f"Worker {worker_rank} terminated. Currently {active_workers} active workers")

    for worker_rank in range(1, size):
        comm.send(None, dest=worker_rank, tag=0)

    end = time.perf_counter()
    print(f"Finished all clustering tasks in {end - start:0.4f} seconds")
    print(f"Total number of clusters: {len(results)}")

    out_dir = os.path.join(out_path, "out")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, "clustering_results.pkl"), 'wb') as f:
        pickle.dump(results, f)

    with open(os.path.join(out_dir,'markers.pkl'), 'wb') as f:
        pickle.dump(markers, f)
    
    # convert the clustering results to a .csv file
    n_cells = sum(len(i) for i in results)
    cl = ['unknown']*n_cells
    for i in range(len(results)):
        for j in results[i]:
            cl[j] = i+1
    res = pd.DataFrame({'cl': cl}, index=adata.obs_names)
    res.to_csv(os.path.join(out_dir,'cl.csv'))
    
    print(f"Manager finished writing clustering results to file")

if __name__ == "__main__":
    if rank == MANAGER_RANK:  # Only the master process executes this
        adata_path = sys.argv[1]
        scvi_path = sys.argv[2]
        out_path = sys.argv[3]
        manager_job_queue(adata_path, scvi_path, out_path)