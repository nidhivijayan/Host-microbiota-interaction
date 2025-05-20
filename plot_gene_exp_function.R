# Define the function
plot_gene_expression <- function(dds_deseq_all, gene_id, metadata, intgroup = c("Treatment"), title = NULL) {
  
  # Step 1: Extract gene counts using plotCounts
  geneCounts <- plotCounts(dds_deseq_all, gene =gene_id, intgroup = intgroup, returnData = TRUE)
  
  # Step 2: Rename the first column to the gene name
  colnames(geneCounts)[1] <- gene_id
  
  # Step 3: Convert rownames to a column and merge with metadata
  geneCountsf <- rownames_to_column(geneCounts, "X")
  geneCountsf <- merge(geneCountsf, metadata, by.x = "X", by.y = "Sample")
  
  # Step 4: Set row names back and remove the "X" column
  geneCountsf <- column_to_rownames(geneCountsf, "X")
  
  # Step 5: Generate the ggplot
  p <- ggplot(geneCountsf, aes(x = Treatment.x, y = get(gene_id), color = Treatment.x)) + #y = get(gene_id)
    scale_y_log10() +
    geom_beeswarm(cex = 3, size = 3) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 60, hjust = 1, size = 12)) +
    scale_color_manual(values = c("#56B4E9", "#b52367","darkgrey"))
  
  # Return the plot
  return(p)
}

library(ggbeeswarm)
library(tibble)
library(ggpubr)
library(ggplot2)
library(DESeq2)
theme_set(theme_bw())
plot_gene_expression(dds_deseq_all, gene_id = "Dhrs3", metadata = metadata, title = "Kmt2e")+theme(axis.text.y=element_text(size=12))+
  stat_compare_means(comparisons = a_my_comparisons, symnum.args = symnum.args, method = "t.test",size=7)+theme_set(theme_bw())+theme(axis.text=element_text(size=12,angle=30,hjust=1))

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

a_my_comparisons <- list(c("HhWT+anti-p40","HhWT+isotype_control_mAb"),c("HhWT+anti-p40","No_treatment_controls"),c("No_treatment_controls","HhWT+isotype_control_mAb"))

c("HhWT+anti-p40","No_treatment_controls")
