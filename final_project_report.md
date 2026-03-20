# EE282 Transcriptomic Analysis of Normal and Glaucomatous Human Optic Nerve Head tissue
Glaucoma is the leading cause of irreversible blindness worldwide, and is characterized by progressive vision loss due to pressure-mediated damage to retinal ganglion cell axons, which all pass through the optic nerve head before innervating visual targets in the brain.
Both neuronal and glial cell populations in this region display profound, often heterogeneous changes in several models of glaucoma.
Therefore, the optic nerve head is an important region of cellular and molecular changes, although the transcriptomic changes in human glaucoma have never been identified at the single cell level.
Here, we performed single cell transcriptomics on optic nerve head tissues from normal and glaucomatous donor eyes.

## Methods
Donor optic nerve heads were collected and frozen less than 10 hours post mortem, and mechanically dissociated to isolate single nuclei. Droplet-based capture and library preparation was performed using the 10X Genomics 3' v4  chemistry.
Illumina sequencing was performed at the UC Irvine core Genomics facility, targetting 300M reads per sample.
Raw sequences were first processed using **CellRanger** and subsequently loaded into **Seurat** for further analysis. The R package **SoupX** was applied to remove ambient RNA and clustered with a resolution of 0.4.
To account for batch effects, the software package **Harmony** was also used.
Clusters were analyzed for top 10 markers and this list was used to assign likely cell type identities to each cluster, based on manual verification.
Astrocytes were further analyzed for DEGs and analyzed for pathway changes using the R package **clusterProfiler**.

## CellRanger
CellRanger is the default software package, provided by 10X Genomics, for analyzing 10X single cell data.
It performs demultiplexing, read alignment, barcode and unique molecular identifier (UMI) processing, feature counting, and generation of a feature barcode matric to be fed into Seurat.

cellranger.sh script:

```
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
                     --create-bam=false \
                     --localcores=4 \
                     --localmem=100
```
### Plot QC metrics
VlnPlot(seu_merged, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "Group", 
        pt.size = 0, 
        ncol = 3)
<img width="661" height="557" alt="QCmetrics" src="https://github.com/user-attachments/assets/9620fef8-36a6-46c0-a011-7a2affc62ba5" />


## Create Seurat object
seurat.sh
```
#!/bin/bash
#SBATCH --job-name=ONH_Mult
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=64G
#SBATCH --partition=standard
#SBATCH --account=gzode_lab
#SBATCH --array=0
#SBATCH --output=seurat_%A_%a.out
#SBATCH --error=seurat_%A_%a.err

#Load libraries, set files and directory

source /dfs6/pub/mnahmou/sc_env/bin/activate

SAMPLES=("ONH_G1a" "ONH_G2a "ONH_N3a" "ONH_N4a")
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

#SoupX configuration
toc <- Read10X(paste0(data_dir, "filtered_feature_bc_matrix/"))
tod <- Read10X(paste0(data_dir, "raw_feature_bc_matrix/"))
sc <- SoupChannel(tod, toc)

#Create temporary clusters for SoupX to work with
temp_seurat <- CreateSeuratObject(counts = toc)
temp_seurat <- NormalizeData(temp_seurat, verbose = FALSE)
temp_seurat <- FindVariableFeatures(temp_seurat, verbose = FALSE)
temp_seurat <- ScaleData(temp_seurat, verbose = FALSE)
temp_seurat <- RunPCA(temp_seurat, verbose = FALSE)
temp_seurat <- FindNeighbors(temp_seurat, dims = 1:20, verbose = FALSE)
temp_seurat <- FindClusters(temp_seurat, resolution = 0.5, verbose = FALSE)
sc <- setClusters(sc, setNames(temp_seurat\$seurat_clusters, colnames(temp_seurat)))

# Remove ambient RNA
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
message("SUCCESS: Analysis complete for ", sample_id)
EOF

## Harmony batch correction
# harmony.sh

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
```

## Cluster cell type annotation
The following R code was applied to the integrated and batch-corrected Seurat object.

```
#Find defining markers and list the top 10 unique markers with highest Log2FC.
cluster_markers <- FindAllMarkers(seu_merged, 
                                  only.pos = TRUE,      
                                  min.pct = 0.25,       
                                  logfc.threshold = 0.25)
top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
top10_genes <- unique(top10_markers$gene)

#Create a dotplot
DotPlot(seu_merged, features = top10_genes) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(title = "Top 10 Marker Genes per Cluster", x = "Genes", y = "Cluster")


seu_merged <- RenameIdents(seu_merged, 
                           "0" = "Oligodendrocytes",
                           "1" = "Astrocytes",
                           "2" = "Photoreceptors",
                           "3" = "Amacrine Cells",
                           "4" = "Low Quality / Dying Cells",
                           "5" = "Reactive Astrocytes",
                           "6" = "ONH Fibroblasts",
                           "7" = "Amacrine Cells",
                           "8" = "Cone Photoreceptors",
                           "9" = "Microglia",
                           "10" = "ON-Bipolar Cells",
                           "11" = "Rod Photoreceptors",
                           "12" = "Retinal Neurons",
                           "13" = "Müller Glia",
                           "14" = "Uncharacterized Glia/Fibroblasts",
                           "15" = "Rod Photoreceptors",
                           "16" = "Cone Photoreceptors",
                           "17" = "Horizontal Cells",
                           "18" = "Amacrine Cells",
                           "19" = "Amacrine Cells",
                           "20" = "OFF-Bipolar Cells",
                           "21" = "OFF-Bipolar Cells",
                           "22" = "Endothelial",
                           "23" = "Retinal Ganglion Cells",
                           "24" = "RPE",
                           "25" = "Amacrine Cells")
DimPlot(seu_merged, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) +
  ggtitle("Annotated ONH Cell Types")
```
<img width="2000" height="750" alt="top10_dotplot" src="https://github.com/user-attachments/assets/62afc2eb-1b35-4784-91e4-3f45c8ead109" />
<img width="1025" height="613" alt="Rplot" src="https://github.com/user-attachments/assets/b382a66b-add7-4767-8204-1c1741a1aec6" />


## DEGs in astrocytes
```
#Create DEG list, sorted by significance
astro_degs <- FindMarkers(seu_merged, 
                          ident.1 = c("ONH_G1a", "ONH_G2a"),
                          ident.2 = c("ONH_N3a", "ONH_N4a"),
                          group.by = "orig.ident",
                          subset.ident = "Astrocytes",
                          min.pct = 0.1, 
                          logfc.threshold = 0.1)
astro_degs <- astro_degs %>% arrange(p_val_adj)
write.csv(astro_degs, file = "~/Desktop/Astrocytes_Glaucoma_vs_Normal.csv", row.names = TRUE)
```

## DEGs in **reactive** astrocytes
```
#Create DEG list, sorted by significance
reactive_astro_degs <- FindMarkers(seu_merged, 
                                   ident.1 = c("ONH_G1a", "ONH_G2a"), 
                                   ident.2 = c("ONH_N3a", "ONH_N4a"), 
                                   group.by = "orig.ident", 
                                   subset.ident = "Reactive Astrocytes", 
                                   min.pct = 0.1, 
                                   logfc.threshold = 0.1)
reactive_astro_degs <- reactive_astro_degs %>% arrange(p_val_adj)
write.csv(reactive_astro_degs, file = "~/Desktop/Reactive_Astrocytes_Glaucoma_vs_Normal.csv", row.names = TRUE)
```

## Pathway analysis on astrocytes for UPregulated genes
```
#Load Packages and DEG list
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

astro_degs <- read.csv("~/Desktop/Astrocytes_Glaucoma_vs_Normal.csv")
colnames(astro_degs)[1] <- "gene"
sig_up_genes <- astro_degs %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% 
  pull(gene)
go_results_up <- enrichGO(gene          = sig_up_genes,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
dotplot(go_results_up, showCategory = 10) + 
  ggtitle("Upregulated Pathways: Glaucoma vs Normal Astrocytes")
```
<img width="765" height="452" alt="astrocyte_upegulated" src="https://github.com/user-attachments/assets/12f5494f-464c-4c90-b572-1f072eef1b3a" />

## Pathway analysis on astrocytes for DOWNregulated genes
```
#Load Packages and DEG list
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

astro_degs <- read.csv("~/Desktop/Astrocytes_Glaucoma_vs_Normal.csv")
colnames(astro_degs)[1] <- "gene"
sig_down_genes <- astro_degs %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>% 
  pull(gene)
go_results_down <- enrichGO(gene          = sig_down_genes,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
dotplot(go_results_down, showCategory = 10) + 
  ggtitle("Downregulated Pathways: Glaucoma vs Normal Astrocytes")
```
<img width="751" height="613" alt="astrocyte_downregulated" src="https://github.com/user-attachments/assets/2e907696-fa29-4372-a71a-9a84eb00b04a" />

## Pathway analysis on reactive astrocytes for UPregulated genes
```
#Load Packages and DEG list
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

reactive_degs <- read.csv("~/Desktop/Reactive_Astrocytes_Glaucoma_vs_Normal.csv")
sig_up_reactive <- reactive_degs %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% 
  pull(gene)
go_reactive_up <- enrichGO(gene          = sig_up_reactive,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = "SYMBOL",  
                           ont           = "BP",      
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)
dotplot(sig_up_reactive, showCategory = 10) + 
  ggtitle("Upregulated Pathways: Glaucoma vs Normal Reactive Astrocytes")
```
<img width="751" height="613" alt="reactive_astrocytes_upregulated" src="https://github.com/user-attachments/assets/de40b416-2f0a-488c-bb8a-3681ff146f55" />

## Pathway analysis on reactive astrocytes for DOWNregulated genes
```
#Load Packages and DEG list
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

reactive_degs <- read.csv("~/Desktop/Reactive_Astrocytes_Glaucoma_vs_Normal.csv")
sig_down_reactive <- reactive_degs %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>% 
  pull(gene)
go_reactive_down <- enrichGO(gene          = sig_down_reactive,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = "SYMBOL",  
                           ont           = "BP",      
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)
dotplot(sig_down_reactive, showCategory = 10) + 
  ggtitle("Downregulated Pathways: Glaucoma vs Normal Reactive Astrocytes")
```
<img width="751" height="613" alt="reactive_astrocytes_downregulated" src="https://github.com/user-attachments/assets/e79bc019-b25b-40a7-90c4-698e598d52d2" />

















