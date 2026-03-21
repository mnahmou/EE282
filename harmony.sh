#!/bin/bash
#SBATCH --job-name=harmony
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=164G
#SBATCH --partition=standard
#SBATCH --account=gzode_lab
#SBATCH --array=0
#SBATCH --output=harmony_%A_%a.out
#SBATCH --error=harmony_%A_%a.err

library(Seurat)
library(harmony)
library(dplyr)

# Load files and metadata
g1a <- readRDS("ONH_G1a_processed.rds")
g2a <- readRDS("ONH_G2a_processed.rds")
n3a <- readRDS("ONH_N3a_processed.rds")
n4a <- readRDS("ONH_N1a_processed.rds")

g1a$Sample <- "ONH_G1a"
g1a$Group  <- "Group_G"
g2a$Sample <- "ONH_G2a"
g2a$Group  <- "Group_G"
n3a$Sample <- "ONH_N3a"
n3a$Group  <- "Group_N"
n4a$Sample <- "ONH_N4a"
n4a$Group  <- "Group_N"

# Combine samples and remove individuals
seu_merged <- merge(x = g1a, 
                    y = c(g2a, n1a, n3a), 
                    add.cell.ids = c("G1a", "G2a", "N1a", "N3a"), 
                    project = "ONH_Project")
seu_merged[["RNA"]] <- JoinLayers(seu_merged[["RNA"]])

rm(g1a, g2a, n1a, n3a)
gc()

#Assign variables and run Harmony
seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)
seu_merged <- ScaleData(seu_merged)
seu_merged <- RunPCA(seu_merged, npcs = 30)
seu_merged <- RunHarmony(seu_merged, 
                         group.by.vars = "Sample", 
                         reduction = "pca", 
                         assay.use = "RNA", 
                         reduction.save = "harmony")

#Create harmonized Seurat object
seu_merged <- FindNeighbors(seu_merged, reduction = "harmony", dims = 1:30)
seu_merged <- FindClusters(seu_merged, resolution = 0.4)
saveRDS(seu_merged, "ONH_integrated_harmony.rds")

echo "Harmony processing done."
