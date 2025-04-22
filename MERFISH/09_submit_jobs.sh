#!/bin/bash
#SBATCH --job-name=process_batches
#SBATCH --output=slurm-logs/slurm-%x.%j.out
#SBATCH --error=slurm-logs/slurm-%x.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --partition=cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dhenze@stanford.edu

# Load modules from Vizgen.sif container
module load anaconda
conda activate Vizgen_2

# Define the paths
ADATAPATH="ABC_cleaned.h5ad"
BASEPATH="/hpc/projects/group.quake/doug/Shapes_Spatial/"
OUTPUTDIR="shape_outputs/"

cd /hpc/mydata/doug.henze/MERFISH/Baysor

# Create the output directory if it doesn't exist
#mkdir -p $OUTPUTDIR
mkdir -p shape_outputs
# Read batch IDs from the adata file
BATCH_IDS=$(python -c "import scanpy as sc; ad = sc.read_h5ad('${ADATAPATH}'); print(' '.join(ad.obs['batchID'].unique()))")

# Array to hold job IDs
JOB_IDS=()

# Submit a separate job for each batch_id
for BATCH_ID in $BATCH_IDS; do
    JOB_ID=$(sbatch --parsable --export=ALL,ADATAPATH=$ADATAPATH,BASEPATH=$BASEPATH,BATCH_ID=$BATCH_ID,OUTPUTDIR=$OUTPUTDIR process_batch.sh)
    JOB_IDS+=($JOB_ID)
done

# Create a comma-separated list of job IDs
JOB_IDS_STR=$(IFS=,; echo "${JOB_IDS[*]}")

# Submit the concatenation job with dependency on the completion of all batch jobs
sbatch --dependency=afterok:$JOB_IDS_STR concatenate_job.sh $OUTPUTDIR feature_vectors_texture.csv
