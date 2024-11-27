from mpi4py import MPI
import queue
import pickle
import time
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import re
import ast

from transcriptomic_clustering.iterative_clustering import onestep_clust, OnestepKwargs

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MANAGER_RANK = 0

def manager_job_queue(adata_path, latent_path, out_path, clust_kwargs): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()

    def preprocess_dict(dict_str):
    # Remove comments using regex
        dict_str_cleaned = re.sub(r"#.*", "", dict_str)
        return dict_str_cleaned
        
    clust_kwargs = preprocess_dict(clust_kwargs)
    clust_kwargs = ast.literal_eval(clust_kwargs)

    means_vars_kwargs = clust_kwargs['means_vars_kwargs']
    highly_variable_kwargs = clust_kwargs['highly_variable_kwargs']
    pca_kwargs = clust_kwargs['pca_kwargs']
    filter_pcs_kwargs = clust_kwargs['filter_pcs_kwargs']
    filter_known_modes_kwargs = clust_kwargs['filter_known_modes_kwargs']
    latent_kwargs = clust_kwargs['latent_kwargs']
    cluster_louvain_kwargs = clust_kwargs['cluster_louvain_kwargs']
    merge_clusters_kwargs = clust_kwargs['merge_clusters_kwargs']

    clust_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )

    adata = sc.read(adata_path)
    print(f"Finished reading in anndata: {adata}")

    if np.max(adata.X) > 100:
        print(f"Raw count data provided")
        print(f"Normlazing total counts to 1e6...")
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        print(f"Finished normalization. max:{np.max(adata.X)}")
    else:
        print(f"Normalized data provided")

    # determine if latent_path is valid
    if os.path.exists(latent_path):
        latent = pd.read_csv(latent_path, index_col=0)
        latent = latent.loc[adata.obs_names]
        adata.obsm['latent'] = np.asarray(latent)
        clust_kwargs.latent_kwargs['latent_component'] = 'latent'
        print(f"Finished reading in latent space: {latent.shape}")
    else:
        print(f"latent_path is invalid or not provided. Using latent_kwargs['latent_component'] for clustering (None for PCA and a str for obsm key)")
    
    min_samples = 4
    random_seed = 2024
    min_samples = 4
    tmp_dir =  os.path.join(out_path,'tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # let the manager do the first clustering
    clusters, markers = onestep_clust(adata, clust_kwargs, random_seed)
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

                new_task = (onestep_clust, (new_adata_path, clust_kwargs, random_seed, idx))
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

                    new_task = (onestep_clust, (new_adata_path, clust_kwargs, random_seed, idx))
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

    print(f"Finished writing clustering results to {out_dir}")
    print(f"Next step: run a final merge of clusters")

if __name__ == "__main__":
    if rank == MANAGER_RANK:  # Only the master process executes this
        adata_path = sys.argv[1]
        latent_path = sys.argv[2]
        out_path = sys.argv[3]
        clust_kwargs = sys.argv[4]
        manager_job_queue(adata_path, latent_path, out_path, clust_kwargs)