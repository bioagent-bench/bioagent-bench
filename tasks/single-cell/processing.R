# Introduction ----
# this script walks through the quality assessment (QA) and analysis of single cell RNA-seq data
# In the 1st 1/2 of the script, we'll practice some basics using a small (~1000 cell) dataset from human peripheral blood mononuclear cells (PBMCs). This dataset comes from the public datasets on the 10X Genomics website: https://www.10xgenomics.com/resources/datasets
# In the 2nd 1/2 of the script, we'll import two separate Seurat objects generated from the spleen of naive and Toxoplasma gondii infected mice, giving us an opportunity to create and analyze an integrated dataset

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("DropletUtils")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("scater")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("Scran")


# Load packages ----
library(tidyverse)
library(DropletUtils) # robust stats package for removing empty droplets from your data
library(Seurat) # a huge, powerful, and popular library for analyzing single cell genomic data
library(scater) # bioconductor package for quality control and visualization for scRNA-seq data
library(scran) # bioconductor package for low level processing of scRNA-seq data
library(Matrix)
library(scales)
library(rjson)
library(R2HTML)
library(DT)

# Import data into R and filter out empty drops ----
# Begin by setting up a new RProject in the folder where you just processed your scRNA-seq data with Kb
# load raw data matrix using the readMM function from the Matrix package
raw_mtx <- readMM('processed/3_quantification/counts_unfiltered/cellranger/matrix.mtx') 
# load genes
genes <- read.csv('processed/3_quantification/counts_unfiltered/cellranger/genes.tsv', sep = '\t', header = F)
# add ensemble gene_ids to the data matrix as rownames
rownames(raw_mtx) <- genes[,1] 
# add cell barcodes as column names
colnames(raw_mtx) <- read.csv('processed/3_quantification/counts_unfiltered/cellranger/barcodes.tsv', sep = '\t', header = F)[,1] 

# use DropletUtils package to get probability that each barcode is a cell
out <- emptyDrops(raw_mtx) 
# set threshold probability for calling a cell
keep <- out$FDR <= 0.05 
# use threshold to remove empty drops
keep[is.na(keep)] <- FALSE
filt_mtx <- raw_mtx[,keep] 

# write out filtered results
write10xCounts('processed/4_analysis/counts_filtered', gene.symbol = genes[,2], filt_mtx, overwrite=T) 

# Generate QC report ----
# this report will contain some useful metrics as well as the traditional log-transformed UMI rank plot (a.k.a. 'waterfall' or 'knee' plot)
# this plot was first described in the Drop-seq paper: - Macosko et al. 2015, DOI:10.1016/j.cell.2015.05.002
# this plot has two important features that we will try to identify:
# 1. 'knee point' - is the point where the signed curvature is minimized. 
# This corresponds to a transition between a distinct subset of barcodes with large totals and the majority of barcodes with smaller totals
# 2. 'inflection point' - is the point on the curve where the first derivative is minimized. 
# This corresponds to the point past which cells cannot reliably be distinguished from background

# source the R script that contains the bc_rank_plot and print_HTML functions we'll use to produce a QC report
# this script comes from Sarah Ennis's github repo here:  https://github.com/Sarah145/scRNA_pre_process
source('./plot_functions.R') 

# load filtered mtx
filt_mtx <- readMM('processed/3_quantification/counts_filtered/matrix.mtx') 

# load run info from JSON files produced by Kb
kb_stats <- c(fromJSON(file = 'processed/3_quantification/inspect.json'), 
              fromJSON(file = 'processed/3_quantification/run_info.json')) 

# determine chemistry version
tech <- grep('processed/3_quantification/10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) 

# make a nice/simple table that summarizes that stats
seq_stats <- data.frame(stat = c('Sequencing technology', 'Number of reads processed', '% reads pseudoaligned', # get sequencing/alignment stats 
                                 '% reads on whitelist'), 
                        value = prettyNum(c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
                                            round(kb_stats$percentageReadsOnWhitelist,2)), big.mark = ','))

# calculate cell stats and save to df
p_cnts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2) 
med_cnts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated number of cells', '% counts in cells', 
                                  'Median counts per cell', 'Median genes per cell', 'Total genes detected'), 
                         value = prettyNum(c(ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
                                             med_genes_cell, tot_genes_detected), big.mark = ','))

# get rank stats
stats <- barcodeRanks(raw_mtx)

# load raw cells
raw_cells <- read.csv('processed/3_quantification/counts_unfiltered/cellranger/barcodes.tsv', header = F, sep ='\t')[,1] 

# load filtered cells
filt_cells <- read.csv('processed/4_analysis/counts_filtered/barcodes.tsv', header = F, sep ='\t')[,1] 

# create barcode rank plot png
bc_rank_plot(stats = stats, 
             raw_cells = raw_cells, 
             filt_cells = filt_cells, 
             save = 'processed/4_analysis/counts_filtered/barcode_rank.png') 

# output a HTML summary of the run
print_HTML(seq_stats = seq_stats, 
           cell_stats = cell_stats, 
           dir = 'processed/4_analysis/counts_filtered', 
           sample_id = NULL)

# Create Seurat object ----
# we'll use the filtered matrix from our emptyDrops analysis to create the Seurat object
# if you're working with data from CellRanger, you will still use the Read10X function below to read in your filtered feature_bc_matrix file
datadir <- 'processed/4_analysis/counts_filtered'
list.files(datadir)

expression_matrix <- Read10X(
  datadir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# actually creating the Seurat Object
pbmc.1k.seurat <- CreateSeuratObject(counts = expression_matrix, min.cells = 3) 

# Now we normalize the Seurat object we just created and find variable features
pbmc.1k.seurat <- NormalizeData(pbmc.1k.seurat, verbose = FALSE) %>% 
pbmc.1k.seurat <- FindVariableFeatures(pbmc.1k.seurat, verbose = FALSE)

# Let's calculate percent of mitochondrial reads
# NOTE: change 'MT' to 'mt' for mouse
pbmc.1k.seurat[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.1k.seurat, pattern = "^MT-") 
# in the violin plot above, features = genes detected, while counts = total molecules detected
# Make violin plot
VlnPlot(pbmc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1)
# Filter your data
pbmc.1k.seurat <- subset(pbmc.1k.seurat, subset = 
                           nCount_RNA < 20000 & 
                           nCount_RNA > 1000 & 
                           nFeature_RNA > 1000 & 
                           percent.mt < 20)
# NOTE: you need to be careful when setting cut-offs that you're not losing unique cell populations

# another QA plot
FeatureScatter(pbmc.1k.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Potential things to look for in the type of QA plot produced above:
# 1. Data points in the bottom LEFT hand quadrant = low genes and UMIs per cell. May represent poor quality cells.
# 2. Data points in the bottom RIGHT hand quadrant = low genes but high UMIs per cell. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# take a look at the top variable features 
top10 <- head(VariableFeatures(pbmc.1k.seurat), 10)
top10

# Plot UMAP ----
# it is standard practice to apply a linear transformation ('scaling') before PCA. For single cell data this includes:
# 1. Shifting the expression of each gene, so that the mean expression across cells is 0
# 2. Scaling the expression of each gene, so that the variance across cells is 1
# This gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(pbmc.1k.seurat)
pbmc.1k.seurat <- ScaleData(pbmc.1k.seurat, features = all.genes)
pbmc.1k.seurat <- RunPCA(pbmc.1k.seurat, npcs = 40, verbose = FALSE)
# What contributes to the PCA?  Let's take a look using the VizDimLoadings function
VizDimLoadings(pbmc.1k.seurat, dims = 1:2, reduction = "pca")
# How many dimensions should we keep? Usually not all 40 
ElbowPlot(pbmc.1k.seurat) #in this case I wouldn't keep more than 20 for the steps below

pbmc.1k.seurat <- RunUMAP(pbmc.1k.seurat, reduction = "pca", dims = 1:20)
pbmc.1k.seurat <- FindNeighbors(pbmc.1k.seurat, reduction = "pca", dims = 1:20)
pbmc.1k.seurat <- FindClusters(pbmc.1k.seurat, resolution = 0.5)
DimPlot(pbmc.1k.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Find cluster-specific genes ----
# generally speaking there are three main ways you can find cluster-specific marker genes with Seurat
# 1. 'FindMarkers' to compare a select cluster to all other cells not in that cluster
# 2. 'FindAllMarkers' to compare EACH cluster to all other cells
# 3. 'FindConservedMarkers' to identify genes conserved (shared) between two defined clusters

# We'll start with FindMarkers, since it allows you to choose exactly which cluster you'd like to focus on.
# Seurat has implemented a speed-saving measure using presto package. Have them install this package and load in the beginning
library(presto) # install with devtools::install_github("immunogenomics/presto")
cluster1.markers <- FindMarkers(pbmc.1k.seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers$pct.diff <- cluster1.markers$pct.1 - cluster1.markers$pct.2
cluster1.markers.df <- as_tibble(cluster1.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits_cluster1 <- cluster1.markers.df %>% arrange(desc(avg_log2FC))
myTopHits_cluster1 <- dplyr::slice(myTopHits_cluster1, 1:20)

# you can make this list of genes into an interactive table
datatable(myTopHits_cluster1, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Cluster 1 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

# plot genes of interest on UMAP
FeaturePlot(pbmc.1k.seurat, 
            reduction = "umap", 
            features = c("GZMB"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)

# now let's try with FindAllMarkers
pbmc.1k.markers <- FindAllMarkers(pbmc.1k.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# let's take the top 10 marker genes for each cluster and plot as a heatmap
top10 <- pbmc.1k.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc.1k.seurat, features = top10$gene)

# Assigning identity to cell clusters  ----
library(SingleR) #automated cell type annotation ('label transfer') using reference data
library(celldex) #a large collection of reference expression datasets with curated cell type labels for use with SingleR package
library(pheatmap)
library(Azimuth) # providing the code below for help with Azimuth installation - Say "yes" to source installation
# remotes::install_github('satijalab/azimuth', ref = 'master')
# depending on your existing package environment, unfortunately this can sometimes be a headache to install
# If azimuth fails to install, you can try the following
# set a higher timeout option for R as follows:
    # options(timeout = max(1000, getOption("timeout")))
    # this high timeout will let you install the following dependency, if R is complaining that it wants it
    # install.packages("BSgenome.Hsapiens.UCSC.hg38")
# now try to install Azimuth again
    # remotes::install_github('satijalab/azimuth', ref = 'master')
# if running Azimuth fails and complains with the error "object 'CRsparse_colSums' not found", this is a documented error (https://github.com/satijalab/seurat/issues/8202#issue-2047511055)
# you can fix this with by reinstalling TFBStools as follows:
    # BiocManager::install("TFBSTools", type = "source", force = TRUE)

# it can also be useful to turn the Seurat object into a singleCellExperiment object, for better interoperability with other bioconductor tools
# two ways to get singleCellExperiment object
# option 1 - use 'read10xCounts' function from DropletUtils package
pbmc.1k.sce <- read10xCounts(datadir)

# option 2 - Seurat allows you to convert directly
pbmc.1k.sce <- as.SingleCellExperiment(pbmc.1k.seurat)

# the singleCellExperiment data structure is easy to work with
rownames(pbmc.1k.sce)
colnames(pbmc.1k.sce)
reducedDims(pbmc.1k.sce)
assays(pbmc.1k.sce)
my.subset <- pbmc.1k.sce[,c(1,2,8)]
rowData(pbmc.1k.sce)$Symbol <- rownames(pbmc.1k.sce)

# OPTION 1: assign identity to cell clusters using public RNA-seq datasets
# To do this, we'll use singleR and celldex (requires an internet connection to connect to ExperimentHub)
ENCODE.data <- BlueprintEncodeData(ensembl = FALSE) #259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE
HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE) #713 microarray samples from the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
DICE.data <- DatabaseImmuneCellExpressionData(ensembl = FALSE) #1561 bulk RNA-seq samples of sorted immune cell populations
ImmGen.data <- ImmGenData(ensembl = FALSE) # 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen)
Monaco.data <- MonacoImmuneData(ensembl = FALSE) #114 bulk RNA-seq samples of sorted immune cell populations that can be found in GSE107011.
MouseRNAseq.data <- MouseRNAseqData(ensembl = FALSE) #358 bulk RNA-seq samples of sorted cell populations
Hemato.data <- NovershternHematopoieticData(ensembl = FALSE) #211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759

predictions <- SingleR(test=pbmc.1k.sce, assay.type.test=1, 
                       ref=Monaco.data, labels=Monaco.data$label.main)

plotScoreHeatmap(predictions)

#now add back to singleCellExperiment object (or Seurat objects)
pbmc.1k.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(pbmc.1k.sce, colour_by = "SingleR.labels")

# OPTION 2: assign identity to cell clusters using Azimuth
# Azimuth using reference 'atlas' scRNA-seq datasets to ID clusters in a query dataset
# Azimuth lets you choose different levels of specificity for cell type annotation
# The RunAzimuth function can take a Seurat object as input
pbmc.1k.seurat <- RunAzimuth(pbmc.1k.seurat, reference = "pbmcref", 
                   query.modality = "RNA", 
                   umap.name = "azimuth.umap")
# visualize UMAP with Azimuth labels
DimPlot(pbmc.1k.seurat, reduction = "umap", group.by  = "predicted.celltype.l1", label = TRUE) + NoLegend()

# Integrate multiple scRNA-seq datasets ----
# To demonstrate integration, we'll leave behind the PBMC dataset we worked with above
# We'll read in two Seurat objects - one generated from the spleen of a untreated mouse (control), and the second from the spleen of mouse infected with Toxoplasma gondii
load("spleen.naive.seurat")
DimPlot(spleen.naive.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

load("spleen.toxoInfected.seurat")
DimPlot(spleen.toxoInfected.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

# since we are now going to work with multiple samples, we need a study design file with our sample metadata
studyDesign_singleCell <- read_tsv("studyDesign_singleCell.txt")

# extract variables of interest
sampleID <- studyDesign_singleCell$sampleID
treatment <- studyDesign_singleCell$treatment

# annotate your seurat objects with as much or as little metadata as you want!
spleen.naive.seurat$treatment <- treatment[1]
spleen.toxoInfected.seurat$treatment <- treatment[2]

# merge the two seurat objects together as a single object (note: this is NOT integration)
spleen_merged <- merge(spleen.naive.seurat, y = c(spleen.toxoInfected.seurat),
                           add.cell.ids = c("spleen.naive.seurat", "spleen.toxoInfected.seurat"),
                           project = "Combined", merge.data = TRUE)

spleen_merged <- split(spleen_merged, f = spleen_merged$treatment)

# Run the standard workflow for visualization and clustering
spleen_merged <- NormalizeData(spleen_merged)
spleen_merged <- FindVariableFeatures(spleen_merged)
spleen_merged <- ScaleData(spleen_merged)
spleen_merged <- RunPCA(spleen_merged)
spleen_merged <- FindNeighbors(spleen_merged, dims = 1:30, reduction = "pca")
spleen_merged <- FindClusters(spleen_merged, resolution = 2, cluster.name = "unintegrated_clusters")
spleen_merged <- RunUMAP(spleen_merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(spleen_merged, reduction = "umap.unintegrated", group.by = c("treatment", "seurat_clusters"))

# Now we actually integrate!
# Importantly, the layered integration works based on how the samples are split - above in the code we split the samples according to treatment so this code will integrate across treatment levels
spleen_integrated <- IntegrateLayers(object = spleen_merged, 
                                     method = CCAIntegration, 
                                     orig.reduction = "pca", 
                                     new.reduction = "integrated.cca",
                                     verbose = FALSE)

# re-join layers after integration - quick
spleen_integrated[["RNA"]] <- JoinLayers(spleen_integrated[["RNA"]])
# repeat clustering and visualization on the integrated seurat object
spleen_integrated <- FindNeighbors(spleen_integrated, reduction = "integrated.cca", dims = 1:30)
spleen_integrated <- FindClusters(spleen_integrated, resolution = 1)
spleen_integrated <- RunUMAP(spleen_integrated, dims = 1:30, reduction = "integrated.cca")
DimPlot(spleen_integrated, reduction = "umap", group.by = c("treatment", "seurat_clusters"))

# let's see what proportion of our total cells reside in each cluster
prop.table(table(Idents(spleen_integrated)))

# remember, we have metadata in this integrated seurat object, so you can use this to split your UMAP
DimPlot(spleen_integrated, reduction = "umap", 
        #split.by = "treatment", # this facets the plot 
        group.by = "seurat_clusters", # labels the cells with values from your group.by variable
        label = TRUE)

# plot genes of interest on UMAP
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = 'Sdc1',
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

# we can plot more than one gene here
my_fav_genes <- c("Cd4", "Cd8a")
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = my_fav_genes,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

# Leveraging cluster identity in your analysis ----
# now let's rerun our cluster identification using SingleR
spleen_integrated.sce <- as.SingleCellExperiment(spleen_integrated)
predictions <- SingleR(test=spleen_integrated.sce, assay.type.test=1, 
                       ref=MouseRNAseq.data, labels=MouseRNAseq.data$label.main)

# now add back to singleCellExperiment object (or Seurat objects)
spleen_integrated.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(spleen_integrated.sce, colour_by = "SingleR.labels")

spleen_integrated2 <- as.Seurat(spleen_integrated.sce, counts = NULL)
DimPlot(spleen_integrated2, reduction = "UMAP", 
        split.by = "treatment", # this facets the plot 
        group.by = "SingleR.labels", # labels the cells with values from your group.by variable
        label = TRUE)

# Now let's try annotation with Azimuth
spleen_integrated <- RunAzimuth(spleen_integrated, reference = "pbmcref", 
                             query.modality = "RNA", 
                             umap.name = "azimuth.umap")

# visualize UMAP with Azimuth labels
DimPlot(spleen_integrated, reduction = "umap", group.by  = "predicted.celltype.l1", label = TRUE) + NoLegend()

# If we now take all the info we get from SingleR and Azimuth, together with our own manual curation using marker genes, we can arrive at a consensus for what each cluster might be.
# keep in mind that this process takes a lot of manual work, and the more detailed you want to get with your clustering, the more work it will be]
new.cluster.ids <- c("B cells", "RBCs", "B cells", "B cells", "CD8+ T cells", "RBCs", "CD4+ T cells", "CD4+ T cells", "Monocytes/Macrophages", "Granulocytes", "B cells", "Monocytes/Macrophages", "Plasma cells", "Granulocytes", "CD8+ T cells", "Monocytes/Macrophages", "NK cells", "Monocytes/Macrophages", "Dendritic cells", "19", "20", "21", "22", "23", "24") 
Idents(spleen_integrated) <- spleen_integrated$seurat_clusters
names(new.cluster.ids) <- levels(spleen_integrated$seurat_clusters)
spleen_integrated <- RenameIdents(spleen_integrated, new.cluster.ids)
# take a look at what you've done
Idents(spleen_integrated)
# we finish by plotting the UMAP with the new cluster labels 
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        label = TRUE)