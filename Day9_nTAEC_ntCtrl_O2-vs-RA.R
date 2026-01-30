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


# Turn condition into a factor with levels for metadata










