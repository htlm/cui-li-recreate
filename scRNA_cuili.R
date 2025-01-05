# Data from https://www.journal-of-hepatology.eu/article/S0168-8278(23)00190-3/fulltext
# 23 patients, scRNA-seq data of liver cells


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
# using both @ and [[]] operators. Right now it has 3 fields per celL: dataset ID,
# number of UMI reads detected per cell (nCount_RNA), and the number of expressed 
# (detected) genes per same cell (nFeature_RNA).

meta <- srat@meta.data
dim(meta)
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)





# note that the data passed to the assay slot has to be a matrix!
tung <- SingleCellExperiment(
  assays = list(counts = as.matrix(tung_counts)),
  colData = tung_annotation
)

# remove the original tables as we don't need them anymore
rm(tung_counts, tung_annotation)

# print the current sce object
tung

##### Access information about SCE objects
rowData(tung)	#Table of gene metadata. There is none in these files
colData(tung)	#Table of cell metadata. Cell-level grouped metadata
class(colData(tung)) # "DFrame" "S4Vectors". DFrame is a type of data.frame used in bioconductor
assay(tung, "counts")	#The assay named “counts.”
class(assay(tung,"counts")) #"matrix" "array" 
reducedDim(tung, "PCA")	#The reduced dimensionality table named “PCA” -- there is no such table rn

# These commands are equivalent
counts(tung)
assay(tung, "counts")

# Investigating batches
unique(tung$batch) # Unique batches in tung
table(tung$batch) # Unique batches with the number of cells in each batch -- each is 96, probably sets of 96 well plates



##### Modifying sce objects
assay(tung, "logcounts") <- log2(counts(tung) + 1) # Use the conventional name 'logcounts'
#Because we named our assay “logcounts,” and this is one of the conventional assay names, we can use the dedicated function to access it:
# first 10 rows and 4 columns of the logcounts assay
logcounts(tung)[1:10, 1:4] # or: assay(tung, "logcounts")[1:10, 1:5]


# Here is a summary of other ways to modify data in SCE objects:
assay(sce, "name") <- matrix	#Add a new assay matrix. The new matrix has to have matching rownames and colnames to the existing object.
rowData(sce) <- data_frame	#Replace rowData with a new table (or add one if it does not exist).
colData(sce) <- data_frame	#Replace colData with a new table (or add one if it does not exist).
colData(sce)$column_name <- values	#Add a new column to the colData table (or replace it if it already exists).
rowData(sce)$column_name <- values	#Add a new column to the rowData table (or replace it if it already exists).
reducedDim(sce, "name") <- matrix	#Add a new dimensionality reduction matrix. The new matrix has to have matching colnames to the existing object.


##### Matrix Statistics
colMeans(counts(tung))
colData(tung)$mean_counts <- colMeans(counts(tung)) # Adding this metadata to our sce object
colData(tung) # Confirm that the new metadata has been added to the sce object

# Other summary statistics
# row (feature) summaries
rowSums(counts(tung))  # sum
rowMeans(counts(tung)) # mean
rowSds(counts(tung))   # standard deviation
rowVars(counts(tung))  # variance
rowIQRs(counts(tung))  # inter-quartile range
rowMads(counts(tung))  # mean absolute deviation

# column (sample) summaries
colSums(counts(tung))  # sum
colMeans(counts(tung)) # mean
colSds(counts(tung))   # standard deviation
colVars(counts(tung))  # variance
colIQRs(counts(tung))  # inter-quartile range
colMads(counts(tung))  # mean absolute deviation



#### Exercises
# Ex 1: Add a new column to colData named “total_counts” with the sum of counts in each cell.
colData(tung)$total_counts  <- colSums(counts(tung))

# Ex 2: Create a new assay called “cpm” (Counts-Per-Million), which contains the result of dividing the counts matrix by the total counts in millions.
assay(tung, "cpm") <- sweep(counts(tung),2, colData(tung)$total_counts/1e6,'/') # 2 = col-wise

# Ex 3: How can you access this new assay?
assay(tung, "cpm") # or equivalently cpm(tung)


#### Subsetting/filtering SCE objects
# subset by numeric index
tung[1:3, ] # the first 3 genes, keep all cells
tung[, 1:3] # the first 3 cells, keep all genes
tung[1:3, 1:2] # the first 3 genes and first 2 cells

# subset by name
tung[c("ENSG00000069712", "ENSG00000237763"), ]
tung[, c("NA19098.r1.A01", "NA19098.r1.A03")]
tung[c("ENSG00000069712", "ENSG00000237763"), c("NA19098.r1.A01", "NA19098.r1.A03")]

# calculate the mean counts per gene
gene_means <- rowMeans(counts(tung))

# print the first 10 values
gene_means[1:10]
gene_means[1:10] > 0.01
# We can use such a logical vector inside [ to filter our data, which will return only the cases where the value is TRUE:
tung[gene_means > 0.01, ]

# A common use case is to retain cells with a certain number of genes above a certain threshold of expression. 
# For this question, we need to break the problem into parts. First let’s check in our counts matrix, 
# which genes are expressed above a certain threshold:

# counts of at least 1
counts(tung) > 0 #T/F matrix
# Because TRUE/FALSE are encoded as 1/0, we can use colSums() to calculate the total number of genes above this threshold per cell:
# total number of detected genes per cell
total_detected_per_cell <- colSums(counts(tung) > 0) # Gene expr counts where total gene expression above threshold (here 0) per cell
tt <- colSums(counts(tung) == 0)

# print the first 10 values
total_detected_per_cell[1:10]

# Say we want cells with at least 5000 detected genes:
tung[, total_detected_per_cell > 5000]

# Common filters:
colSums(counts(tung)) > x	#Total counts per cell greater than x.
colSums(counts(tung) > x) > y	#Cells with at least y genes having counts greater than x.
rowSums(counts(tung)) > x	#Total counts per gene greater than x.
rowSums(counts(tung) > x) > y	#Genes with at least y cells having counts greater than x.

# Ex 1: Create a new object called tung_filtered which contains:
#       cells with at least 25000 total counts
#       genes that have more than 5 counts in at least half of the cells

cell_filter <- colSums(counts(tung)) >= 25000 
      # not total_detected_per_cell >= 25000 because 
      # total_detected_per_cell <- colSums(counts(tung) > 0)  This tells you how many genes have counts above 0 in each cell, not actual exp lvls
      #             cell_filter <- colSums(counts(tung)) >= 25000 The parentheses are in different places
gene_filter <- rowSums(counts(tung) > 5) > ncol(tung)/2

tung_filtered <- tung[gene_filter, cell_filter]


##### Visual Data Exploration
# General Framework
# ggplot(data = <data.frame>, 
#        mapping = aes(x = <column of data.frame>, y = <column of data.frame>)) +
#   geom_<type of geometry>()

# Create the input dataframe
cell_info <- as.data.frame(colData(tung))

# Plot the counts distribution in each sc-seq batch
ggplot(data = cell_info, aes(x = batch, y = total_counts)) +
  geom_violin(fill = 'brown') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Plot the counts distribution in each sc-seq batch now using ggcells
  # Note that we can now plot directly from the tung sce object
ggcells(tung, aes(x = batch, y = total_counts)) + 
  geom_violin(fill = 'orange') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Plotting the gene expression for a single gene across each batch
ggcells(tung, aes(x = batch, y = ENSG00000198938), exprs_values = "logcounts") + 
  geom_violin(fill = 'coral2') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Ex 1: Make a scatterplot showing the relationship between the mean and the 
#       variance of the raw counts per cell. (Bonus: also colour the cells by batch.) 
tung$sds_counts <- colSds(counts(tung))
colData(tung)$var_counts <- colVars(counts(tung))
ggcells(tung, aes(x = mean_counts, y = sds_counts, label = batch), show.legend = TRUE) + 
  geom_point(aes(colour = batch)) +
  theme_bw()
ggcells(tung, aes(x = mean_counts, y = var_counts, label = batch), show.legend = TRUE) + 
  geom_point(aes(colour = batch)) + 
  theme_bw()


##### CellRanger Aside
library(DropletUtils)

# importing the raw count data
sce <- read10xCounts("data/pbmc_1k_raw")

# importing the pre-filtered count data
sce_filtered <- read10xCounts("data/pbmc_1k_filtered")

