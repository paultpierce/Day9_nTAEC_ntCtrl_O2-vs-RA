## Load DESeq2 library

library(DESeq2)
library(tidyverse)


# Load in data and metadata of sample info

counts <- read.delim(file = "Comp8_Day9_ntCtrl_O2-vs-RA_RawGeneCounts.csv", sep = ",")
head(counts)

metadata <- read.delim(file = "design.csv", sep = ",")
head(metadata)

# Ensure that rownames have GeneIds as name for DESeq analysis

head(counts)
rownames(counts) = counts$GeneIds
head(counts)

# Remove unused GeneIds and Gene columns from counts

genes = counts[ , c("GeneIds", "Gene")]
counts = counts[ , -c(1, 2)]
head(counts)


# Ensure that rownames have sample IDs  as name in metadata matrix for DESeq analysis

head(metadata)
rownames(metadata) = metadata$sample_ID
head(metadata)

# Check to see if rownames from metadata and colnames of counts match for creating DESeq object 

all(rownames(metadata) == colnames(counts))

# Turn condition into a factor with levels for metadata

metadata$condition = factor(metadata$condition, levels = c("O2", "RA"))
metadata$condition


# Spot checking expression for oxidative stress gene (HMOX1)

gene_id = genes$GeneIds[genes$Gene == "SFTPC"]
gene_counts = counts[gene_id, ]

gene_data = cbind(metadata, counts = as.numeric(gene_counts))


ggplot(gene_data, aes(x = condition, y = counts, fill = condition)) +
        geom_boxplot() +
        theme_bw(base_size = 14) +
        xlab(NULL) +
        labs(title = paste0(gene_id, " raw counts by condition"))


# Creating dds object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)











