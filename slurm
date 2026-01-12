#!/bin/bash -l
# Full-node Snakemake job on Raven

#SBATCH -J pre_process         # job name
#SBATCH -o logs/snakemake_%j.out     # stdout
#SBATCH -e logs/snakemake_%j.err     # stderr
#SBATCH -D /u/yacheng/projects/clear_bd  # working directory

#SBATCH --nodes=1                    # one node
#SBATCH --ntasks-per-node=1          # one task on the node
#SBATCH --cpus-per-task=72           # use all 72 cores on the node
#SBATCH --mem=200000                 # 200 GB RAM
#SBATCH --time=01:00:00              # walltime: 1 hours

module purge

# Load whatever you need for conda (adjust if different on Raven)
# Example for user-local Miniconda:
source ~/miniconda3/etc/profile.d/conda.sh
source ~/.bashrc
# Activate the environment that has snakemake installed
conda activate /u/yacheng/conda-envs/snakemake

mkdir -p logs

# Run Snakemake using all cores on this node
snakemake \
  --snakefile Snakefile \
  --configfile snakemake/config.yaml \
  -q \
  --cores 72 \
  --keep-going

#run second time to make sure all dependencies are updated
snakemake \
  --snakefile Snakefile \
  --configfile snakemake/config.yaml \
  -q \
  --cores 72 \
  --keep-going