#!/bin/bash
#SBATCH --job-name=concatenate_vectors
#SBATCH --output=slurm-logs/slurm-%x.%j.out
#SBATCH --error=slurm-logs/slurm-%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dhenze@stanford.edu

# Load modules from Vizgen.sif container
module load anaconda
conda activate Vizgen_2

# Run the Python script for concatenating feature vectors
python concatenate_feature_vectors.py $1 $2

