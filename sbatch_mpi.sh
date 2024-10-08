#!/bin/bash
#SBATCH --job-name=mpi
#SBATCH --partition=celltypes
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --time=100:00:00
#SBATCH --mem=200G
#SBATCH --mail-user=dan.yuan@alleninstitute.org
#SBATCH --mail-type=END,FAIL

# Things to edit in this .sh file: 
# 1. Number of nodes and size to use. for over 1 million cells, I recommend using 11 nodes with each 500G memory. It will take 2 hours to finish. 
#  - 5 200G nodes are enough for 20k cells
#  - If HPC is busy, reduce the number of nodes. If requested n nodes, change the number of nodes for the workers to use in the last line of code to (n-1)
# 2. the absolute path to the h5ad file, scvi latent space (.csv file), and a output path where the results will be saved. The .h5ad file should have raw counts or normalized counts in adata.X. The algorithm will automatically detect if the counts are raw or normalized. If it is raw counts, it will do the normalization. This will be noted in the log file
# 3. the absolute path to the manager script and the worker script

adata_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/HMBA_analysis/iscANVI_mapping/troubleshoot/adata_query_MHGlut.h5ad" # norm or counts in X is fine.
scvi_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/HMBA_analysis/iscANVI_mapping/troubleshoot/run3/scvi_09062024/scvi_latent_integrated.csv" 
out_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc" # output directory

manager_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/iterative_clustering_mpi_manager.py"
worker_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/iterative_clustering_mpi_worker.py"

source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/rapids_singlecell

module load mpi/mpich-3.2-x86_64

# Navigate to the output directory, where the log and out files will be saved
# Check if the out directory exists. if not, create it
if [ ! -d "$out_dir" ]; then
    # If the directory doesn't exist, create it
    mkdir -p "$out_dir"
    echo "Directory '$out_dir' created."
fi
cd "$out_dir"

export PYTHONPATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/tc_latent:$PYTHONPATH

time mpiexec -n 1 sh -c "python \"$manager_script\" \"$adata_dir\" \"$scvi_dir\" \"$out_dir\"> manager_output.log 2> manager_error.log" : -n 4 sh -c "python \"$worker_script\" > worker_output.log 2> worker_error.log" # somehow worker_output.log is trucated, but its ok.
