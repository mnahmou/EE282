#!/bin/bash
#SBATCH --job-name=cell_ranger
#SBATCH --account=gzode_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=164G
#SBATCH --time=96:00:00
#SBATCH --partition=standard
#SBATCH --output=cellranger_multi_%j.out
#SBATCH --error=cellranger_multi_%j.err

cd /pub/mnahmou/sc/
module load cellranger/8.0.1
SAMPLES=("ONH_N1a" "ONH_G1a" "ONH_G2a" "ONH_N4a")
for SAMPLE in "${SAMPLES[@]}"
do
    cellranger count --id="${SAMPLE}" \
                     --transcriptome="${PWD}/refdata-gex-GRCh38-2024-A" \
                     --fastqs="${PWD}/input" \
                     --sample="${SAMPLE}" \
                     --create-bam=true \
                     --localcores=8 \
                     --localmem=124
