#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
conda activate ee282
# Run busco on fly genome
busco -i hifi_fly_assembly.fasta -o busco_fly_eval -m genome -l diptera_odb10 -c 16
