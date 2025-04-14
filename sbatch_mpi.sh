#!/bin/bash
#SBATCH --job-name=mpi
#SBATCH --partition=celltypes
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --mail-user=dan.yuan@alleninstitute.org
#SBATCH --mail-type=END,FAIL

# Things to edit in this .sh file: 
# 1. Number of nodes and size to use. for over 1 million cells, I recommend using 11 nodes with each 500G memory. It will take 2 hours to finish. 
#  - 5 200G nodes are enough for 20k cells
#  - If HPC is busy, reduce the number of nodes. If requested n nodes, change the number of nodes for the workers to use in the last line of code to (n-1)
# 2. the absolute path to the h5ad file, scvi latent space (.csv file), and a output path where the results will be saved. The .h5ad file should have raw counts or normalized counts in adata.X. The algorithm will automatically detect if the counts are raw or normalized. If it is raw counts, it will do the normalization. This will be noted in the log file
# 3. the absolute path to the manager script and the worker script
# 4. the clustering parameters. The default parameters are set to match the bigcat default parameters. You can change them as needed.

adata_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/test/subsample_data/test_subsample_1000k_count.h5ad"
latent_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/test/subsample_data/scvi_latent_1000k/scvi_latent_1millionCells.csv"
out_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/archive_versions/v3/test_1M_2"

manager_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/hicatMPI/iterative_clustering_mpi_manager.py"
worker_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/hicatMPI/iterative_clustering_mpi_worker.py"

# The following are the default parameters for the clustering. You can change them as needed.
clust_kwargs="{
    'means_vars_kwargs': {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, change to 0.6931472 to match to the bigcat default
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    },
    'highly_variable_kwargs': {
        'max_genes': 4000 # originally 3000, change to 4000 to match to the bigcat default
    },
    'pca_kwargs': {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    },
    'filter_pcs_kwargs': {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    },
    # if not using known_modes, set filter_known_modes_kwargs to None or an empty dict. Only applies for PCA
    'filter_known_modes_kwargs': {
        'known_modes': 'log2ngene', 
        'similarity_threshold': 0.7
    },
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    'latent_kwargs': {
        'latent_component': None # None (default to run PCA) or obsm key such as "X_scVI"
    },
    'cluster_louvain_kwargs': {
        'k': 15, # number of nn, originally 150, change to 15 to match to the bigcat default
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus, originally 8
        'resolution': 1.0 # resolution of louvain for taynaud method
    },
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, change to 10 to match to the bigcat default
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 1, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, change to 100 to match to the bigcat default
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, change to 0.6931472 to match to the bigcat default
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, change to 4 to match to the bigcat default
        'de_method': 'ebayes'
    }
}"

source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc

if [ ! -d "$out_dir" ]; then
    # If the directory doesn't exist, create it
    mkdir -p "$out_dir"
    echo "Directory '$out_dir' created."
fi
cd "$out_dir" # navigate to the output directory, where the log and out files will be saved

export PYTHONPATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering:$PYTHONPATH

time mpiexec -n 1 sh -c "python \"$manager_script\" \"$adata_dir\" \"$latent_dir\" \"$out_dir\" \"$clust_kwargs\"> manager_output.log 2> manager_error.log" : -n 9 sh -c "python \"$worker_script\" \"$out_dir\" 2> worker_error.log" 