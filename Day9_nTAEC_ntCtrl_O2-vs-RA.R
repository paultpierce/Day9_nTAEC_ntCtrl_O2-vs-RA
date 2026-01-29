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
count_data <- read_csv("Comp8_Day9_ntCtrl_O2-vs-RA_RawGeneCounts.csv")

# Create sample metadata (sample info) including tech replicates
sample_info <- read_csv("design.csv")

# 