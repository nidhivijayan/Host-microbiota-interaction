# https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0710-2#Sec2
#vsd_all from Fig1
keep_all <- rowSums(counts(dds_deseq_all) >= 20)>=3
dds_filtered <- dds_deseq_all[keep_all, ]
vsd_all <- vst(dds_filtered)

res_4_cor1<-results(dds_filtered, contrast = c("Treatment", "HhWT+isotype_control_mAb", "No_treatment_controls"),alpha=0.05)
res_4_cor2<-results(dds_filtered, contrast = c("Treatment", "HhWT+isotype_control_mAb", "HhWT+anti-p40"),alpha=0.05)
res_4_cor3<-results(dds_filtered, contrast = c("Treatment", "No_treatment_controls", "HhWT+anti-p40"),alpha=0.05)

padj_cutoff <- 0.05
logFC_cutoff <- 2 

# Identify significantly UPregulated genes in each treatment
up_genes_Hh <- rownames(res_4_cor1[which(res_4_cor1$padj < padj_cutoff & res_4_cor1$log2FoldChange > logFC_cutoff), ])
up_genes_Hh2 <- rownames(res_4_cor2[which(res_4_cor2$padj < padj_cutoff & res_4_cor2$log2FoldChange > logFC_cutoff), ])
up_C<-rownames(res_4_cor3[which(res_4_cor3$padj < padj_cutoff & res_4_cor3$log2FoldChange > logFC_cutoff), ])

# Identify genes that are UP in Control (i.e., Control > Hh, Control > anti-p40)
control_vs_Hh <- rownames(res_4_cor1[which(res_4_cor1$padj < padj_cutoff & res_4_cor1$log2FoldChange < -logFC_cutoff), ])  # Control > Hh
up_antip40_1 <- rownames(res_4_cor2[which(res_4_cor2$padj < padj_cutoff & res_4_cor2$log2FoldChange < -logFC_cutoff), ])  # Control > anti-p40
up_antip40_2 <- rownames(res_4_cor3[which(res_4_cor3$padj < padj_cutoff & res_4_cor3$log2FoldChange < -logFC_cutoff), ])  # Control > anti-p40

# Merge Control-upregulated genes (allowing overlap with anti-p40)
up_genes_control <- unique(c(control_vs_Hh, up_C)) 

up_antiP40<-unique(c(up_antip40_1, up_antip40_2)) 

all_genes2 <- unique(c(up_antiP40, up_genes_control, up_genes_Hh ,up_genes_Hh2))
up_control<-as.data.frame(up_genes_control)
View(up_control)
write.table(up_control, "up_control.txt")
up_antiP40_df<-as.data.frame(up_antiP40)
write.table(up_antiP40_df, "up_anti-p40.txt")
up_Hh<-as.data.frame(up_genes_Hh)
write.table(up_Hh, "up_Hh.txt")

View(all_genes2)
all_genes22<-as.data.frame(all_genes2)
up_genes_C_df<-as.data.frame(up_genes_C)
up_genes_C_df<-column_to_rownames(up_genes_C_df,"up_genes_C")

all_genes22<-column_to_rownames(all_genes22,"all_genes2") #This fives 1072 genes
genes_4_corr <- merge(as.data.frame(all_genes22), as.data.frame(counts(dds_filtered,normalized =TRUE)), by = 'row.names', sort = FALSE)
up_genes_control_df<-merge(as.data.frame(up_genes_C_df), as.data.frame(counts(dds_filtered,normalized =TRUE)), by = 'row.names', sort = FALSE)
View(up_genes_control_df)
genes_4_corr<-column_to_rownames(genes_4_corr,"Row.names")

genes_4_corr[]<-lapply(genes_4_corr, as.integer)

genes_4_corr<- read.delim("genes_4_corr.txt",check.names = FALSE,row.names ="ID")


dds_corr<-DESeqDataSetFromMatrix(countData =genes_4_corr,colData = metadata, design=~Treatment)
dds_deseq_corr <- DESeq(dds_corr)
vsd_corr <- rlog(dds_deseq_corr)

rna_counts <- assay(vsd_corr)  # Extract matrix
View(rna_counts)

otu_table_log<-read.delim("OTU_table_log_summed_by_family.txt",check.names = FALSE,row.names = "Family")
View(otu_table_log)
common_samples <- intersect(colnames(otu_table_log), colnames(rna_counts))
otu_table_log <- otu_table_log[, common_samples]
rna_counts <- rna_counts[, common_samples]
write.table(rna_counts,"rna_counts_lfc2.txt")
rna_counts <- read.delim("rna_counts_ed.txt",check.names = FALSE,row.names = "ID") #840 genes after removing lnc genes
library(psych)
# Run Spearman correlation for all pairs
cor_results3 <- corr.test(t(otu_table_log), t(rna_counts), method="spearman", adjust="fdr", ci=FALSE)
cor_matrix <- cor(t(otu_table_log), t(rna_counts), method="spearman") #To run without p value
cor_results <- rcorr(t(otu_table_log),t(rna_counts), type="spearman")

View(cor_matrix)
# Extract correlation coefficients and p-values
# cor_matrix <- cor_results2$r
# pval_matrix <- cor_results2$p

cor_matrix <- cor_results3$r
pval_matrix <- cor_results3$p
cor_matrix2<-read.delim("cor_matrix_3.txt",row.names = "X",check.names = FALSE)
pval_matrix2<- read.delim("pval_matrix_3.txt",row.names = "X",check.names=FALSE)
p_threshold <- 0.01
cor_threshold <- 0.70

# Create a logical matrix where TRUE means the correlation meets the threshold
significant_mask2 <- abs(cor_matrix2) >= cor_threshold & pval_matrix2 <= p_threshold

# Identify rows and columns with at least one significant correlation (TRUE in mask)
keep_rows <- rowSums(significant_mask2) > 0
keep_cols <- colSums(significant_mask2) > 0

# Subset the original correlation matrix to keep only rows/columns with significant correlations
cor_matrix_filtered2 <- cor_matrix2[keep_rows, keep_cols]

# View the cleaned correlation matrix
View(cor_matrix_filtered2) #132 ASVs and 709 genes

# Create a logical matrix where TRUE means the correlation is greater than abs(0.7)
above_threshold <- abs(cor_matrix_filtered2) > cor_threshold

# Count the number of values in each column that exceed the threshold
col_counts <- colSums(above_threshold)

# Keep rows where at least 3 values exceed the threshold
cor_matrix_filtered_final_columns <- cor_matrix_filtered2[, col_counts >= 4] #this was 3 : for 0.8 corr: 17 genes, 63 ASVs

# View the filtered correlation matrix
View(cor_matrix_filtered_final_columns) #has 545 genes, and 132 ASVs #for 0.85 corr: 7 genes, 31 ASVs

cor_matrix_filtered_final_columns<- as.matrix(cor_matrix_filtered_final_columns)
#go to Fig4_p2 to filter rows
# Cluster rows
row_hclust <- hclust(dist(cor_matrix_filtered_final_columns), method = "ward.D2")
row_order <- row_hclust$order  # Get order

# Cluster columns
col_hclust <- hclust(dist(t(cor_matrix_filtered_final_columns)), method = "ward.D2")
col_order <- col_hclust$order  # Get order

cor_matrix_reordered <- cor_matrix_filtered_final_columns[row_order, col_order]

# Plot reduced correlation matrix
library(corrplot)
corrplot(cor_matrix_reordered, addrect = 3, method = "circle", is.corr = FALSE)

corrplot(cor_matrix_filtered_final_columns, addrect = 3, method='shade',order='AOE',diag=FALSE, is.corr = FALSE)

corrplot(cor_matrix_reordered, method = "color", order = "hclust", hclust.method = "ward.D2")

cor_edges <- melt(cor_matrix_filtered_final_columns)
# Filter out weak correlations (adjust threshold as needed)
cor_edges <- cor_edges[abs(cor_edges$value) > 0.7, ]

# Create graph from edge list
network <- graph_from_data_frame(cor_edges, directed = FALSE)

# Plot network
plot(network, edge.width = abs(E(network)$value) * 3, vertex.size = 5)

#Export to cytoscape ==
# Create a node data frame
nodes <- data.frame(
  id = unique(c(edge_list$Family, edge_list$Gene)),  # All unique nodes
  type = ifelse(unique(c(edge_list$Family, edge_list$Gene)) %in% edge_list$Gene, "Gene", "Family")  # Label node type
)

# Save as CSV
write.csv(nodes, "nodes_family_top.csv", row.names=FALSE)

# Create edge table
edges <- edge_list %>%
  rename(source = Family, target = Gene)  # Rename for Cytoscape compatibility

# Save as CSV
write.csv(edges, "edges_top.csv", row.names=FALSE)


library(pheatmap)
pheatmap(cor_matrix, color = colorRampPalette(c("blue", "white", "red"))(50))

plot(otu_table_log["Bacteroides", ], rna_counts["Reg3g", ], 
     xlab="Bacteroides (log-abundance)", 
     ylab="Reg3g Expression", 
     main="Spearman Correlation")

threshold <- 0.8
strong_cor <- cor_matrix * (abs(cor_matrix) > threshold)  # Set weak correlations to 0

library(reshape2)

# Convert matrix to long format
cor_long <- melt(cor_matrix)

# Filter for strong correlations (excluding self-correlations r = 1)
strong_cor_df <- cor_long[abs(cor_long$value) > threshold, ]
View(strong_cor_df)

N <- 20  # Number of strongest correlations to keep
top_cor <- head(strong_cor_df[order(-abs(strong_cor_df$value)), ], N)
top_cor$value <- as.numeric(top_cor$value)

ggplot(data = top_cor, aes(x=Var2, y=Var1,
                                   fill=value))+  geom_tile()
#Correlate gene to ASV==
genes_of_interest <- c("Itln1", "Reg3g")

# Extract correlation values for the selected genes
cor_selected <- cor_matrix[, top_cor]

View(cor_matrix)

# Filter for strong correlations (absolute correlation > 0.6) for each gene
strong_with_GeneX <- apply(cor_selected, 2, function(x) x[abs(x) > 0.6])

# Print results
strong_with_GeneX

# strong_with_GeneX <- cor_matrix[, "Itln1"][abs(cor_matrix[, "Itln1"]) > 0.6]

View(strong_with_GeneX)

# Convert list to data frame, filling missing values with NA
library(dplyr)
heatmap_data <- bind_rows(strong_with_GeneX, .id = "Target_Gene") %>%
  column_to_rownames(var = "Target_Gene")

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix[is.na(heatmap_matrix)] <- 0 
# Check structure
str(heatmap_matrix)

library(pheatmap)

pheatmap(heatmap_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Strong Correlations with Itln1 and Reg3g")


#Correlate ASV to gene==
bacteria_of_interest <- c("Campilobacterota;Campylobacteria;Campylobacterales;Helicobacteraceae;Helicobacter")

# Extract correlation values for the selected genes
cor_selected_bac <- cor_matrix[bacteria_of_interest,]

# Filter for strong correlations (absolute correlation > 0.6) for each gene
strong_with_Hel<- cor_matrix["Campilobacterota;Campylobacteria;Campylobacterales;Helicobacteraceae;Helicobacter",][abs(cor_matrix["Campilobacterota;Campylobacteria;Campylobacterales;Helicobacteraceae;Helicobacter",]) > 0.8]


# strong_with_GeneX <- cor_matrix[, "Itln1"][abs(cor_matrix[, "Itln1"]) > 0.6]

View(strong_with_GeneX)

# Convert list to data frame, filling missing values with NA
library(dplyr)
heatmap_data_hel <- bind_rows(strong_with_Hel, .id = "Bacteria") %>%
  column_to_rownames(var = "Bacteria")

# Convert to matrix
heatmap_matrix_hel <- as.matrix(strong_with_Hel)
heatmap_matrix_hel[is.na(heatmap_matrix_hel)] <- 0 
# Check structure
str(heatmap_matrix)

if (ncol(heatmap_matrix) == 1) {
  heatmap_matrix <- t(heatmap_matrix)  # Transpose it
}

pheatmap(heatmap_matrix_hel, 
         color = colorRampPalette(c("blue", "white", "red"))(50),cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Strong Correlations with Helicobacter")



# Convert to matrix for heatmap
gene_subset <- matrix(strong_with_GeneX, nrow = 1, dimnames = list("Itln1", names(strong_with_GeneX)))

# Create heatmap
pheatmap(gene_subset, color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = TRUE, cluster_rows = FALSE, 
         main = "Strong Correlations with Cdc45")



# Remove self-correlations (diagonal) and weak correlations
cor_matrix_filtered <- cor_matrix
cor_matrix_filtered[abs(cor_matrix_filtered) < 0.8] <- NA  # Set weak correlations to NA

# Convert long format to wide format
cor_matrix <- reshape2::acast(top_cor, Var1 ~ Var2, value.var = "value")
cor_matrix[is.na(cor_matrix)] <- 0  # Replace NA with 0

# Plot heatmap
pheatmap(cor_matrix, color = colorRampPalette(c("blue", "white", "red"))(50),
         na_col = "grey", cluster_rows = TRUE, cluster_cols = TRUE)


library(DESeq2)

