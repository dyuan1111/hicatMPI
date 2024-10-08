from mpi4py import MPI
import anndata as ad
import os

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

MANAGER_RANK = 0

def worker_job():
    while True:
        # Receive job from manager
        task = comm.recv(source=MANAGER_RANK, tag=MPI.ANY_TAG)
    
        if task is None:
            # print(f"Worker {rank} received termination signal.")
            # comm.send(None, dest=MANAGER_RANK, tag=0)
            break  # No more jobs, exit loop

        func, args = task
        adata_path, kwargs, random_seed, idx = args

        # Simulate processing (replace with actual computation)
        adata = ad.read_h5ad(adata_path)
        n = adata.shape[0]
        clusters, markers = func(adata, kwargs, random_seed)

        # convert the indices back to the original indices
        original_clusters = []
        for cluster in clusters:
            original_clusters.append([idx[i] for i in cluster])

        # Signal completion to master and request new task
        comm.send((original_clusters, markers, n), dest=MANAGER_RANK) 
        # clusters is a list of lists of cell indices
        # markers is a set of gene names

        # # Print result
        # sizes = [len(cluster) for cluster in clusters]
        # print(f"Worker {rank} finished clustering {adata.shape[0]} cells into {len(clusters)} clusters of sizes {sizes}")

        # delete the temporary anndata file
        os.remove(adata_path)

if __name__ == "__main__":
    worker_job()
