install.packages("WGCNA")
library(WGCNA)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)

counts2<-read.table("merged_counts.txt", header = TRUE,sep = "\t")
colnames(counts2) <- ifelse(grepl("^X", colnames(counts2)), gsub("^X", "", colnames(counts2)), colnames(counts2))
colnames(counts2) <- gsub("\\.", "-", colnames(counts2))
col<-c("GeneID", "4890-MB-1", "4890-MB-2", "4890-MB-3", "4890-MB-4", "4890-MB-5", "4890-MB-6", "4890-MB-7", "4890-MB-8", "4890-MB-9", "4890-MB-10", "4890-MB-11", "4890-MB-12", "4890-MB-13", "4890-MB-14", "4890-MB-15", "4890-MB-16", "4890-MB-17", "4890-MB-18", "4890-MB-19", "4890-MB-20", "4890-MB-21", "4890-MB-22", "4890-MB-23", "4890-MB-24", "4890-MB-25", "4890-MB-26", "4890-MB-27", "4890-MB-28", "4890-MB-29", "4890-MB-30", "4890-MB-31", "4890-MB-32", "4890-MB-33", "4890-MB-34", "4890-MB-35", "4890-MB-36")
counts2<-counts2[,col]
View(counts2)
col_sel = names(counts2)[-1]     # Get all but first column name
mdata <- counts2 %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )

# Optional step ---  order the groups in the plot.
# mdata$group = factor(mdata$group,
#                     levels = c("B", "B_L1", ....))  #<= fill the rest of this in

# ==== Plot groups (Sample Groups vs RNA Seq Counts) to identify outliers # Not needed
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)

wpn_vsd <- getVarianceStabilizedData(dds_deseq_all)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#>  0.00000  0.00000  0.00000  0.08044  0.03322 11.14529

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

input_mat = t(expr_normalized)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#From another turotial ==
adjacency = adjacency(data, power=6, type="unsigned")

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
View(module_df)

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$SampleID = row.names(MEs0) #Change treatment to SampleID

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-Treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=Treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr") #Black and brown asso with p40 and +blue also the nno treatment controls, turqiose and green were associated with inflammation

# pick out a few modules of interest here
modules_of_interest = c("black", "brown", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")+scale_color_manual(values = c("black", "brown", "turquoise"))

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")

edgelist <- read.table("edgelist.tsv", header = FALSE)
g <- graph_from_data_frame(edgelist, directed = FALSE)
threshold <- 0.7
g_subset <- delete_edges(g, E(g)[correlation < threshold])
summary(g_subset)

# Calculate degree centrality
degree_centrality <- degree(g_subset)

# Sort and select top hub genes (e.g., top 10)
top_hub_genes <- names(sort(degree_centrality, decreasing = TRUE))[1:10]

# View top hub genes
top_hub_genes

library(clusterProfiler)

# Assume you have a list of hub genes
hub_genes <- top_hub_genes
View(hub_genes)
hub_genes_symbol <- hub_genes
go_results <- enrichGO(gene = hub_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

# Map symbols to Entrez IDs
gene_mapping_biomart <- getBM(attributes = c('ensembl_gene_id', 'interpro_short_description'),
                              values = hub_genes_symbol,  # This is your vector of Ensembl IDs
                              mart = ensembl)
# GO Enrichment Analysis (for human genes, use the appropriate OrgDb)
View(gene_mapping_biomart)

# View the results
summary(go_results)


# Remove all rows with less than n counts across all samples, where n=#samples
low_count_mask <- rowSums(counts) < ncol(counts)

sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))
# add a colorbar along the heatmap with sample condition
num_conditions <- nlevels(metadata$Colitis)
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
metadata$Colitis <- as.factor(metadata$Colitis)
cond_colors <- pal[as.integer(metadata$Colitis)]

log_counts <- log2(counts + 1)

x = melt(as.matrix(log_counts))

colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + geom_density()
heatmap.2(cor(log_counts),RowSideColors=cond_colors,
          trace='none', main='Sample correlations (log2-transformed)')
# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:
# mod <- model.matrix(~0+samples$condition+samples$batch)
mod <- model.matrix(~0+metadata$Treatment)

# make model terms easier to work with
colnames(mod) <- levels(metadata$Treatment)

fit <- lmFit(log_counts, design=mod)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(samples$condition), 2))                                                                                                                               

comparisons <- list()                                                                                                                                          
for (i in 1:nrow(condition_pairs)) {                                                                                                                                     
  comparisons[[i]] <- as.character(condition_pairs[i,])                                                                                                      
}    

# vector to store de genes
sig_genes <- c()

# iterate over the contrasts, and perform a differential expression test for
# each pair
for (conds in comparisons) {
  # generate string contrast formula, "infLM24 - infLM4"
  contrast_formula <- paste(conds, collapse=' - ')
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
  contrast_fit <- contrasts.fit(fit, contrast_mat)
  eb <- eBayes(contrast_fit)
  
  # Grab highly ranked genes; this is a pretty stringent p-value cutoff, but
  # it serves to limit the total number of genes we will use for this
  # tutorial
  sig_genes <- union(sig_genes, 
                     rownames(topTable(eb, number=Inf, p.value=0.005)))
}

# Filter out genes which were not differentially expressed for any contrast
log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]
