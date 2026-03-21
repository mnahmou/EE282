#!/bin/bash
#SBATCH --job-name=ONH_Mult
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=standard
#SBATCH --account=gzode_lab
#SBATCH --array=0-3
#SBATCH --output=seurat_%A_%a.out
#SBATCH --error=seurat_%A_%a.err

#Load libraries, set files and directory
source /dfs6/pub/mnahmou/sc_env/bin/activate
SAMPLES=("ONH_G1a" "ONH_G2a" "ONH_N3a" "ONH_N4a")
CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
/dfs6/pub/mnahmou/sc_env/bin/Rscript - <<EOF

library(Seurat)
library(SoupX)
library(conflicted)
library(dplyr)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

sample_id <- "$CURRENT_SAMPLE"
data_dir <- paste0(sample_id, "/outs/")

#SoupX to remove ambient RNA
toc <- Read10X(paste0(data_dir, "filtered_feature_bc_matrix/"))
tod <- Read10X(paste0(data_dir, "raw_feature_bc_matrix/"))
sc <- SoupChannel(tod, toc)
sc <- autoEstCont(sc)
adj_counts <- adjustCounts(sc)

# Create and filter the seurat object
onh <- CreateSeuratObject(counts = adj_counts, project = sample_id)
onh[["percent.mt"]] <- PercentageFeatureSet(onh, pattern = "^MT-")
onh <- subset(onh, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

#Process and save 
onh <- NormalizeData(onh)
onh <- FindVariableFeatures(onh)
onh <- ScaleData(onh)
onh <- RunPCA(onh)
onh <- RunUMAP(onh, dims = 1:20)
saveRDS(onh, paste0(sample_id, "_processed.rds"))

EOF
