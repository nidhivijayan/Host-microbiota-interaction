counts<-read.table("merged_counts_v3.txt", header = TRUE,row.names = "GeneID",check.names = FALSE)
col<-c("4890-MB-1", "4890-MB-2", "4890-MB-3", "4890-MB-4", "4890-MB-5", "4890-MB-6", "4890-MB-7", "4890-MB-8", "4890-MB-9", "4890-MB-10", "4890-MB-11", "4890-MB-12", "4890-MB-13", "4890-MB-14", "4890-MB-15", "4890-MB-16", "4890-MB-17", "4890-MB-18")

counts<-counts[,col]

total_reads <- colSums(counts)  # Adjust if needed, skipping annotation columns

df <- data.frame(Sample = names(total_reads), Reads = total_reads)
df <- merge(df, metadata, by = "Sample")

library(ggplot2)
df$Sample <- factor(df$Sample, levels = c("4890-MB-1", "4890-MB-2", "4890-MB-3", "4890-MB-4", "4890-MB-17", "4890-MB-18", "4890-MB-5", "4890-MB-6", "4890-MB-7", "4890-MB-8", "4890-MB-9", "4890-MB-10", "4890-MB-11", "4890-MB-12", "4890-MB-13", "4890-MB-14", "4890-MB-15", "4890-MB-16"))  # Customize order

reads_plot<-ggplot(df, aes(x = Sample, y = Reads, fill = Treatment)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Total Reads per Sample", x = "Sample", y = "Total Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("#02adba","#b52367","grey" ))

reads_plot

# which genes are driving the PCA plot
pca_res <- prcomp(t(assay(vsd_all)))  # t() transposes the data, rows = samples, columns = genes
pca_loadings <- pca_res$rotation 

# For the first principal component (PC1 pr PC2)
pc1_loadings <- pca_loadings[, 1] 
pc2_loadings <- pca_loadings[, 2]
top_genes_pc2 <- order(abs(pc2_loadings), decreasing = TRUE)[1:10]  # Top 10 genes

# Plot the top contributing genes' expression
top_gene_names <- rownames(vsd_all)[top_genes_pc2]
pheatmap(assay(vsd_all)[top_genes_pc1,], trace="none", dendrogram="column", col=brewer.pal(9, "Blues"))

pc2_df <- data.frame(
  gene = rownames(vsd_all)[top_genes_pc2],
  score = pc2_loadings[top_genes_pc2]
)

ggplot(pc2_df, aes(x = reorder(gene, score), y = score)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  labs(x = "Gene", y = "PC1 Score") +theme(axis.text.x=element_text(size=14))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +theme_minimal()
