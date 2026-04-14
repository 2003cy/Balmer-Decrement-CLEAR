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
#SBATCH --time=00:30:00              # walltime: 1 hours

module purge

mkdir -p logs

# Run Snakemake using all cores on this node
pixi run snakemake \
  --snakefile Snakefile \
  --configfile snakemake/config.yaml \
  -q \
  --cores 72 \
  --keep-going

#run second time to make sure all dependencies are updated
pixi run snakemake \
  --snakefile Snakefile \
  --configfile snakemake/config.yaml \
  -q \
  --cores 72 \
  --keep-going

#collect results from Snakemake
pixi run python 7_collect_results.py