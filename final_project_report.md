# EE282 Transcriptomic Analysis of Normal and Glaucomatous Human Optic Nerve Head tissue
## Introduction
Glaucoma is the leading cause of irreversible blindness worldwide, and is characterized by progressive vision loss due to pressure-mediated damage to retinal ganglion cell axons (Ehrlrich et al 2022), all of which pass through the optic nerve head before innervating visual targets in the brain.
Both neuronal and glial cell populations in this region display profound, often heterogeneous changes in several models of glaucoma (Mazumder et al 2023, Cameron et al 2024).
Therefore, the optic nerve head is an important region of cellular and molecular changes, although the transcriptomic changes in human glaucoma in the optic nerve head have never been profiled using single cell/single nuclei RNA sequencing.
Here, we performed single nuclei transcriptomics and pathway analysis on optic nerve head tissues from normal and glaucomatous donor eyes. Two samples from normal donor eyes (referred to as "N3a" and "N4a") and two samples from glaucomatous donor eyes (referred to as "G1a" and "G2a") were considered for this analysis.

Astrocytes, the most abundant nonneuronal cell in the nervous system are important mediators of neuronal disease and become reactive in every neurodegenerative disease, displaying hallmarks such as proliferation, migration, hypertrophy, and changes in gene expression (Cameron et al 2024). Astrocyte reactivity, or astrogliosis is a heterogeneous, context-dependent response to cellular damage, with some reactive astrocytes actively driving neuronal loss through active complement signalling (Liddelow et al 2017). Therefore we focus on changes in astrocytes for this report, aiming to understand their cellular processes that are changing in glaucomatous tissues.

## Methods
Donor optic nerve heads were collected and frozen less than 10 hours post-mortem, and mechanically dissociated to isolate single nuclei. Droplet-based capture and library preparation was performed using the 10X Genomics 3' v4  chemistry.
Illumina sequencing was performed at the UC Irvine Genomics High-Throughput Facility, targetting 300M reads per sample.
Raw sequences were first processed using **CellRanger** and subsequently loaded into **Seurat** for further analysis. The R package **SoupX** was applied to remove ambient RNA and clustered with a resolution of 0.4.
To account for batch effects, the software package **Harmony** was also used.
Clusters were analyzed for top 10 markers and this list was used to assign likely cell type identities to each cluster, based on manual verification.
Astrocytes were further analyzed for DEGs and analyzed for pathway changes using the R package **clusterProfiler**.

## Results
### CellRanger
CellRanger is the default software package, provided by 10X Genomics, for analyzing 10X single cell data (Zheng, et al 2017). The specific package 'count' was implemented for these data.
It performs demultiplexing, read alignment, barcode and unique molecular identifier (UMI) processing, feature counting, and generation of a feature barcode matrix to be fed into downstream analysis modules such as Seurat.
FastQC was proposed in the analysis proposal but was not implemented here because it is not optimized for single cell/single nuclei data. Applying FastQC may produce false sample failures due to the unique architecture of single cell libraries, especially the presence of nucleotide barcodes.

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
                     --create-bam=true \
                     --localcores=8 \
                     --localmem=124
```

The CellRanger output websummary identified quality issues with both normal samples: "Low fraction reads in cells. Ideal > 70%. Application performance may be affected. Many of the reads were not assigned to cell-associated barcodes. This could be caused by high levels of ambient RNA or by a significant population of cells with a low RNA content..."

#### Plot QC metrics
VlnPlot(seu_merged, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "Group", 
        pt.size = 0, 
        ncol = 3)
<img width="661" height="557" alt="QCmetrics" src="https://github.com/user-attachments/assets/9620fef8-36a6-46c0-a011-7a2affc62ba5" />

Violin plots show the distribution of three important metrics: feature(cell) count, transcript count, and the percetange reads of mitochondrial transcripts.
Normal samples have significantly fewer cells and transcripts while having higher mitochondrial content, all indicative of poor sample quality.
SoupX and Harmony will be run to account for the observed batch effects in the following analysis using Seurat in R.

Next, the Seurat package was used to cluster the cellranger feature barcode matrix, (Zheng, et al 2017), implementing SoupX (Young et al, 2020) to clean the data by removing ambient RNA content. When cells are dissociated, some RNAs are floating in the 'soup' and are captured in droplets along with nuclei. SoupX estimates and removes these background RNAs so they do not influence downstream analysis.
Cells were filtered to contain less than 5% mitochondrial DNA and only keep cells with over 200, but not more than 5000 genes per cell, avoiding empty droplets and doublets (droplets with two or more cells) in the Seurat object creation.
Importantly, the deplyr package called in this script has conflicting command names with Seurat so the package 'conflicted' was loaded to manage "filter" and "select" commands. 

### Create Seurat object
seurat.sh
```
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
```
### Harmony batch correction 
Harmony "projects cells into a common embedding" to make samples more comparable (Korsunsky et al 2019). Samples were assigned into biological groups to preserve biological differences while correcting for sample to sample variability.

The following script for harmony.sh was run:

```
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

### Cluster cell type annotation
The following R code was applied to the integrated and batch-corrected Seurat object, using the 'FindAllMarkers' function, which compares gene expression in all clusters combined to every individual cluster and was set to only consider positively regulated genes, increased in expression by 25% in at least 25% of the cells in each cluster. The list of Top10 genes was pasted into Google Gemini in order to assign the most likely cell identity. AI-driven annotations were spot checked for accuracy, especially for "astrocyte" and "reactive astrocyte" clusters by plotting canonical marker genes over the UMAP to confirm colocalization with identified clusters. All expected cell types are present, with astrocytes, reactive astrocytes, and oligodendrocytes constituting half of all cells (50.5%) (24742 cells out of 48949 total cells).

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


Astrocytes are abundant cells in this tissue and are particularly interesting across models of neurodegenerative disease because of their contributions to maintaining neuronal physiology as well as their dynamic responses to various injury cues in their conversion to reactive states. Therefore, astrocyte clusters presented here "Astrocytes" and "Reactive Astrocytes" were the focus of differentially expressed genes (DEG) and pathway analysis performed using the package clusterProfiler which references the Gene Ontology (GO) database of pathway definitions (Yu et al 2012). 

DEGs were computed by the FindMarkers function, selecting for a log fold 2 change of at least 0.1 and in at least 10% of cells.

### DEGs in astrocytes
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

### DEGs in **reactive** astrocytes
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

DEG lists were filtered to select only cells with greater than 25% change in either direction , with a p-value of 0.05 or less. Gene Ontology (GO) enrichment plots are shown below:

### Pathway analysis on astrocytes for UPregulated genes
```
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

### Pathway analysis on astrocytes for DOWNregulated genes
```
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

### Pathway analysis on reactive astrocytes for UPregulated genes
```
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
dotplot(sig_down_reactive, showCategory = 10) + ggtitle("Downregulated Pathways: Glaucoma vs Normal Reactive Astrocytes")
```
<img width="727" height="557" alt="downregulated_reactive astrocytes" src="https://github.com/user-attachments/assets/f4ff4b29-3d5e-4729-8c88-c2b81d8dae39" />

## Findings and discussion
In comparison to control, glaucomatous astrocytes upregulate protein localization to organelles, Wnt signalling, and response to oxidative stress and downregulate RNA splicing, microtubule-based movement, and cilium assembly.
GO terms upregulated in reactive astrocytes include Wnt signaling, actin filament organization, and macroautophagy, while downregulated pathways include RNA splicing and protein transport.

Astrocyte reactivity is itself a heterogeneous process, for which future analysis can address by subclustering all astrocyte populations and identifying pathway changes in each subcluster. Previous reports have identified C3 positive and proliferative astrocytes as critical influences of neuronal survival after injury or in disease (Cameron et al 2024), which we expect to see in this dataset as well. In addition, identification of new markers for different reactive subtypes may be useful for the understanding of reactive astrocyte biology as a whole.

While the data obtained from donor tissues produces cell clusters and identifies associated changes in gene expression, the inherent nature of human sample variability, especially time before freezing, may influence the reliablity of conclusions in this study. Analysis of more samples will be necessary (N = 6-8) in order to identify robust changes in these tissues and motivate follow up mechanistic studies.

Additionally, the set of molecular changes in other cell clusters of interest including oligodendrocytes and Muller glia may be explored in subsequent analyses.

## References

Cameron EG, Nahmou M, Toth AB, Heo L, Tanasa B, Dalal R, Yan W, Nallagatla P, Xia X, Hay S, et al. 2024. A molecular switch for neuroprotective astrocyte reactivity. Nature 626: 574–582.

Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh P, Raychaudhuri S. 2019. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16: 1289–1296.

Liddelow SA, Guttenplan KA, Clarke LE, Bennett FC, Bohlen CJ, Schirmer L, Bennett ML, Münch AE, Chung W-S, Peterson TC, et al. 2017. Neurotoxic reactive astrocytes are induced by activated microglia. Nature 541: 481–487.

Mazumder AG, Julé AM, Sun D. 2023. Astrocytes of the optic nerve exhibit a region-specific and temporally distinct response to elevated intraocular pressure. Mol Neurodegener 18: 68.

Satija R, Farrell JA, Gennert D, Schier AF, Regev A. 2015. Spatial reconstruction of single-cell gene expression data. Nat Biotechnol 33: 495–502.

Tham Y-C, Li X, Wong TY, Quigley HA, Aung T, Cheng C-Y. 2014. Global prevalence of glaucoma and projections of glaucoma burden through 2040: a systematic review and meta-analysis. Ophthalmology 121: 2081–2090.

Young MD, Behjati S. 2020. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience 9: giaa151.

Yu G, Wang L-G, Han Y, He Q-Y. 2012. clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. OMICS 16: 284–287.

Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, et al. 2017. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8: 14049.














