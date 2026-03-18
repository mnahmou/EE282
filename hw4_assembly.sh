#!/usr/bin/env bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=standard

#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err

source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
conda activate ee282

# Run hifiasm
hifiasm \
-o hifi_fly_assembly \
-t 16 \
ISO_HiFi_Shukla2025.fasta.gz
