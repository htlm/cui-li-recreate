# Data from https://www.journal-of-hepatology.eu/article/S0168-8278(23)00190-3/fulltext
# 23 patients, scRNA-seq data of liver cells
# The following code is adapted from https://www.singlecellcourse.org/introduction-to-single-cell-rna-seq.html
# to suit Cui and Li et al, Single-cell atlas of the liver myeloid compartment before 
# and after cure of chronic viral hepatitis.


# load the libraries
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(Seurat)

#load data
expression_matrix <- ReadMtx(
  mtx = "counts/matrix.mtx.gz", 
  features = "counts/features.tsv.gz", feature.column = 1,
  cells = "counts/barcodes.tsv.gz"
)

srat <- CreateSeuratObject(counts = expression_matrix)
srat
str(srat)

# Meta.data is the most important field for next steps. It can be accessed 
# using both @ and [[]] operators. Right now it has 3 fields per cell: dataset ID,
# number of UMI reads detected per cell (nCount_RNA), and the number of expressed 
# (detected) genes per same cell (nFeature_RNA).

meta <- srat@meta.data
dim(meta)
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

# Looking for mitochondrial and ribosomal genes
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
head(srat)

# Looking at the data
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

# nFeature_RNA = number of genes detected in each cell 
# nCount_RNA = total number of molecules detected within a cell
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nFeature_RNA", feature2 = "percent.rb")

# QC from the supplemental material: 
# - less than 25% mitochondrial genes
# - between 500 and 7500 genes per cell
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 , 'Low_nFeature', 'Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA > 7500, 'High_nFeature', srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 25,
                       paste('High_percent.mt',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])

# Normalization from supp material: Normalize so gene expr units in each cell sum to 10k
# normalizedata normalizes to 10k by default
srat <- NormalizeData(srat)

# Get an overview of highly variably expressed genes:
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
most.var.10 <- head(VariableFeatures(srat), 10)
most.var.plot.2k <- VariableFeaturePlot(srat)
LabelPoints(plot = most.var.plot.2k, points = most.var.10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave('plots/cuili_most_var_2k.png', plot = most.var.plot.2k)

# Regnerating Fig S1B - Normalized gene expression of housekeeping genes
# Genes of interest: ACTB, B2M, GAPDH
srat[["old.ident"]] <- Idents(object = srat)
srat <- StashIdent(object = srat, save.name = "old.ident")
srat <- SetIdent(srat, value='cuili')

housekeeping.plot <- VlnPlot(srat, features = c("ACTB","B2M","GAPDH"),pt.size = 0.1, same.y.lims=TRUE) & 
  theme(plot.title = element_text(size=10))
ggsave('plots/cuili_housekeeping.png', plot = housekeeping.plot)

# From supp: removed cells enriched for lymphocyte markers to focus on the myeloid compartment
# lymphocyte markers: CD3, CD4, CD8, CD19 (Bcells), CD16 (Nk cells) Ref: https://www.uncmedicalcenter.org/
gene.names <- sort(c(rownames(srat)))
cd.gene.names <- gene.names[grep('CD', gene.names)]
# CD3, CD8, CD16 all don't seem to be present in the downloaded dataset
lymphocyte.plot <- VlnPlot(srat, features = c("CD4","CD19"),pt.size = 0.1, same.y.lims=TRUE) & 
  theme(plot.title = element_text(size=10))
ggsave('plots/cuili_lymphocyte.png', plot = lymphocyte.plot)


