# Extract normalized counts
norm_counts <- counts(dds_filtered, normalized = TRUE)

# Define target gene and specific genes
target_gene <- "Reg4"  # Example target gene
specific_genes <- c("Cdx2", "Sox9", "Gata6","Gli1","Atf2")  # Genes to check for coexpression

# Subset only these genes
expr_data <- norm_counts[c(target_gene, specific_genes), ]

cor_matrix <- cor(t(expr_data), method = "spearman")  # Use "pearson" or "spearman"
print(cor_matrix)

library(ggplot2)

# Example: Scatter plot for Ifng vs GeneA
five<-ggplot(as.data.frame(t(expr_data)), aes(x = Reg4, y = Atf2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +  # Add correlation stat
  theme_minimal() +
  ggtitle("Coexpression: REG4 vs ATF2")

first
sec
three
four
five

ggarrange(first, sec, three, four, five, labels="A","B","C","D","E")
