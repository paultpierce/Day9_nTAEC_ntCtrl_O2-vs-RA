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


# Add collapseReplicates condition to combine tech replicates

ddsColl <- collapseReplicates(dds, groupby = dds$bio_ID, renameCols = FALSE)
dds <- ddsColl


# Pre-filter genes based on number of counts

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]


# Determine size factors to use for normalization during DESeq

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Extract normalized counts

normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)


# Combine gene names to normalized_counts dataframe and write.csv out for AG

normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df <- merge(genes, normalized_counts_df, by = "row.names")
normalized_counts_df$Row.names <- NULL

write.csv(normalized_counts_df, file = "Day9_ntCtrl_O2-vs-RA_NormData_Paul.csv", row.names = FALSE)


# Quality assessment of normalized count data with heatmaps

library(ggrepel)
library(pheatmap)
library(RColorBrewer)


vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
View(vsd_cor)

pheatmap(vsd_cor, annotation = select(metadata, condition), 
         main = "Hierarchical heatmap analysis by condition")



# QA of normalized counts via PCA

plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA by condition") +
  geom_text_repel(aes(label = sample_ID))



# DESeq analysis of dds object

dds <- DESeq(dds)


# DESeq model - dispersion testing

plotDispEsts(dds)


# Compare expression by conditions O2 vs RA

res <- results(dds, contrast = c("condition", "O2", "RA"))
head(res)


# Convert res DESeq into data frame and combine with genes to get Gene names

res_df <- as.data.frame(res)
head(res_df)
head(genes)

res_df = merge(genes, res_df, by = "row.names")
res_df$Row.names <- NULL
head(res_df)


# Order DESeq results by p adj value

res_df_ordered_padj <- res_df[order(res_df$padj), ]
head(res_df_ordered_padj)


# Cont. filtering but by raw pval value and log2FC
filtered_data_DEGs_pval <- res_df %>% 
  filter(res_df$pvalue < 0.1)

filtered_data_DEGs_pval <- filtered_data_DEGs_pval %>% 
  filter(abs(filtered_data_DEGs_pval$log2FoldChange) > 0.5)


# Cont. filtering but by padj value and log2FC
filtered_data_DEGs_padj <- res_df %>% 
  filter(res_df$padj < 0.1)

filtered_data_DEGs_padj <- filtered_data_DEGs_padj %>% 
  filter(abs(filtered_data_DEGs_padj$log2FoldChange) > 0.5)



# Visualizations

library(EnhancedVolcano)

EnhancedVolcano(res_df, lab = rownames(res_df), 
                x = "log2FoldChange", y = "pvalue")









