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
pThr <- 0.05 # adj. p-val threshold
logFCThr <- 1 # log2 fold-change 
baseMeanThr <- 5 # filtering expression level of raw counts
cpmThr <- 1 # copy-per-million


# Set factors for condition
sample_info$condition <- factor(sample_info$condition)
class(sample_info$condition)

# Had to eliminate gene IDs due to mismatch in ncol in count_data != nrow in sample_info
count_data <- count_data[ , -1]

# DESeq object creation and importing of raw count_data with sample_info
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

######## ddsColl <- collapseReplicates(dds, groupby = dds$bio_ID)
######## Ask Abhrajit about the above due to perceived technical replicates in raw counts file


# Setup factor levels for condition in relation to dds object
dds$condition <- factor(dds$condition, levels = c("O2", "RA"))


# Filter genes based on number of counts from value = baseMeanThr
keep <- rowSums(counts(dds)) >= baseMeanThr
dds <- dds[keep, ]

# Perform stat test/analysis to ID DEGs
dds <- DESeq(dds)

deseq_result <- results(dds)
head(deseq_result)

# Plot MA for DESeq result before transforming into data frame
plotMA(deseq_result, ylim = c(-4, 4))



# Change DESeq object to R object as a data frame
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)


# Order DESeq results by p value
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue), ]
head(deseq_result_ordered)



# Cont. filtering but by p adj value and log 2 FC as set by values pThr & logFCThr, respectively
filtered_data_DEGs <- deseq_result %>% 
                        filter(deseq_result$pvalue < pThr)

filtered_data_DEGs <- filtered_data_DEGs %>% 
                        filter(abs(filtered_data_DEGs$log2FoldChange) > logFCThr)





