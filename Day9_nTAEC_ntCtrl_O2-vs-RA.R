### Load packages needed

# Differential expression analysis
library(DESeq2)

# Data handling & plotting
library(tidyverse)

# Log fold change shrinkage
library(apeglm)

# Annotation database human
library(org.Hs.eg.db)

# Visualizations
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)


# Import count data (raw counts file)
count_data <- read.csv("Comp8_Day9_ntCtrl_O2-vs-RA_RawGeneCounts.csv", header = TRUE, row.names = 1)
colnames(count_data)
head(count_data)

# Create sample metadata (sample info) including tech replicates
sample_info <- read.csv("design.csv", header = TRUE, row.names = 1)
colnames(sample_info)
head(sample_info)

# Variables for downstream DEG selection criteria
pThr <- 0.05 # for adj. p-val threshold
logFCThr <- 1 # for log2 fold-change 
baseMeanThr <- 10 # for filtering expression level of raw counts
cpmThr <- 1 # for copy-per-million


# Set factors for condition
sample_info$condition <- factor(sample_info$condition)
class(sample_info$condition)

# Had to eliminate gene IDs due to mismatch in ncol in count_data != nrow in sample_info
count_data <- count_data[ , -1]

# DESeq object creation and importing of raw count_data with sample_info
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)


# Setup factor levels for condition in relation to dds object
dds$condition <- factor(dds$condition, levels = c("O2", "RA"))












