#!/bin/bash
#SBATCH --job-name=process_batch
#SBATCH --output=slurm-logs/slurm-%x.%j.out
#SBATCH --error=slurm-logs/slurm-%x.%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --partition=cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dhenze@stanford.edu

# Load modules from Vizgen.sif container
module load anaconda
conda activate Vizgen_2

# Run the Python script for processing batches
python process_images.py $ADATAPATH $BASEPATH $BATCH_ID $OUTPUTDIR
