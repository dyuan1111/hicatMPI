# Accelerating Transcriptomic Clustering with Distributed Computing

This project implements a **dynamic, asynchronous clustering algorithm** using the **Message-Passing Interface (MPI)** to distribute clustering tasks across multiple HPC nodes. The method efficiently launches clustering jobs in parallel, assigning pending clustering tasks from the queue as soon as a HPC node completes its current task. This improves upon the original recursive clustering algorithm implemented in the Python package transcriptomic_clustering (see Background below for more information), significantly reducing the time to cluster large datasets. 

### Background
The transcriptomic clustering Python package that uses a **scVI latent space** can be found here: [transcriptomic_clustering](https://github.com/AllenInstitute/transcriptomic_clustering/tree/hmba/tc_latent). It is the Python version of the R package [scrattch.hicat](https://github.com/AllenInstitute/scrattch.hicat), both of which perform clustering recursively (depth-first search). The recursive approach can take significant time for large datasets. For instance, clustering 1 million cells can take ~2 days.

### Distributed Clustering Model
This project replaces the depth-first search (DFS) recursive method with a **dynamic, asynchronous manager-worker model**:
- **Manager Node**: Oversees the clustering process, distributing jobs to worker nodes and managing the queue of tasks.
- **Worker Nodes**: Independently perform clustering tasks and, once a job finishes, immediately taking a job from the queue and performing clustering again.
- **Asynchronous Task Distribution**: The system doesn’t wait for all clustering jobs at the same hierarchy to finish before moving on. Instead, as soon as a worker node completes its current job, it directly moves to the next pending task in the queue.

## How It Works

![Process Illustration](images/mpiTC.jpeg)

- **Step 1**: The manager node performs the initial clustering task. Once that finishes, the manager node appends all subsequent clustering jobs for each cluster into the job queue. For each available worker node, the manager assigns a job from the queue.
- **Step 2**: Each worker node performs its assigned clustering task.
- **Step 3**: Once a worker finishes, it sends the result to the manager and immediately takes a pending job from the queue. The manager evaluate the clusters:
  - for clusters that cannot be further clustered (clustered into 1 cluster), the results are added to the final results.
  - for clusters that can be furtehr clustered (clustered into > 1 clusters), the subclusters are added to the queue for further clustering.
- **Step 4**: The process continues until the job queue is empty and all nodes have terminated.

## Benchmarking Results

The table below shows the comparison of run time between the original recursive clustering method and the dynamic, asynchronous MPI-based method.

| Number of Cells | Recursive        | MPI (This Method)      |
|-----------------|------------------|-------------------------|
| 10,000          | 12 minutes       | 9 minutes               |
| 80,000          | 63 minutes       | 18 minutes              |
| 1,000,000       | 47 hours 41 minutes | 2 hours 15 minutes    |

### Notes on Node Utilization:
- Clustering 1.2 million cells were achieved using 10 nodes with 500G memory each. It took 2-3 hours to finish. 

## Prerequisites

Before running the project, ensure you have the following installed:
- Python 3
- MPI from [mpi4py](https://mpi4py.readthedocs.io/en/stable/mpi4py.html)
- Required Python package:
  - `transcriptomic_clustering`

---

## Running the Code

1. Modify the `sbatch_mpi.sh` script to specify the number of nodes and any required configuration for your HPC environment.
2. Submit the job to your HPC using:
   ```bash
   sbatch sbatch_mpi.sh

