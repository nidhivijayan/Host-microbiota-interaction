
# ===================== Load Data =====================
#  Inputting all samples
# ==========================================
counts<-read.table("merged_counts_v2.txt", header = TRUE,row.names = "GeneID",check.names = FALSE)
metadata<-read.delim("rnaseq_metadata_v2.txt",header=TRUE)
View(metadata)
View(counts)

#order the rows with SampleID in metadata the same as columns in counts matrix
#col<-c("4890-MB-1", "4890-MB-2", "4890-MB-3", "4890-MB-4", "4890-MB-5", "4890-MB-6", "4890-MB-7", "4890-MB-8", "4890-MB-9", "4890-MB-10", "4890-MB-11", "4890-MB-12", "4890-MB-13", "4890-MB-14", "4890-MB-15", "4890-MB-16", "4890-MB-17", "4890-MB-18", "4890-MB-19", "4890-MB-20", "4890-MB-21", "4890-MB-22", "4890-MB-23", "4890-MB-24", "4890-MB-25", "4890-MB-26", "4890-MB-27", "4890-MB-28", "4890-MB-29", "4890-MB-30", "4890-MB-31", "4890-MB-32", "4890-MB-33", "4890-MB-34", "4890-MB-35", "4890-MB-36")
col<-c("4890-MB-1", "4890-MB-2", "4890-MB-3", "4890-MB-4", "4890-MB-5", "4890-MB-6", "4890-MB-7", "4890-MB-8", "4890-MB-9", "4890-MB-10", "4890-MB-11", "4890-MB-12", "4890-MB-13", "4890-MB-14", "4890-MB-15", "4890-MB-16", "4890-MB-17", "4890-MB-18")

counts<-counts[,col]
dds<-DESeqDataSetFromMatrix(countData = counts,colData = metadata, design=~Treatment)
dds_deseq_all <- DESeq(dds)
keep_all <- rowSums(counts(dds_deseq_all) >= 1)>=2 #removes counts with > 1 in less than 2 samples. You can play around with this. 
dds_filtered <- dds_deseq_all[keep_all, ]
rld_all<-rlog(dds_filtered) #shrink
#OR
vsd_all <- vst(dds_filtered) 

pcaData<-plotPCA(vsd_all,intgroup=c("Treatment"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

group_sample<-metadata$Treatment
group_sample<-factor(group_sample)

pca_all <- ggplot(pcaData, aes(PC1, PC2, color = group_sample, fill = group_sample)) +
  scale_fill_manual(values = c("#02adba","#b52367","grey" )) +
  scale_color_manual(values = c("#02adba","#b52367","grey")) +  # Ensures colors match
  geom_point(size = 4, stroke = 0.6, shape = 21, color = "black") +  # Ensure shape supports fill
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title=element_text(size=13),
        legend.text = element_text(size = 12), 
        legend.key.size = unit(2, "lines"))  

pca_all #All above can go to github

#
sampleDists <- dist(t(assay(vsd_all)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(assay(vsd_all))

colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

#rownames(sampleDistMatrix) <- paste(vsd_all$Treatment, vsd_all$p40)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df <- as.data.frame(colData(dds_deseq_all)[, c("Treatment")])

column_to_rownames(metadata,"Sample")
#colnames(df2) <- "Treatment"  # Set correct column name for annotation
rownames(df2) <- colnames(assay(vsd_all))  # Set row names to match sample IDs in vsd_all

annotation_colors <- list(
  Treatment = c(
    "HhWT+anti-p40" = "#56B4E9",
    "HhWT+anti-p40+VPI" = "#E69F00",
    "HhWT+isotype_control_mAb"="#F0E442",
    "HhWT+isotype_control_mAb+VPI"="#009E73",
    "No_treatment_controls"="grey"),
  anti.p40=c(
    "No"="green","Yes"="magenta"),Colitis=c("No"="black","Yes"="pink"),
  Cdifficile=c("No"="black","Yes"="pink"),Sex=c("Male"="black","Female"="pink"))

p1<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, annotation_col = df,annotation_colors = annotation_colors)


# ===================== Add Ensembl =====================
##  load biomart ==========================================
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

View(ensembl)
# Convert final results .csv file into .txt file
results_csv <- "Hh2-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Hh2-replaceoutliers-results.txt"
a <- read.table(results_txt, head=TRUE)
View(a)

b2_part1 <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                  values=a$X,
                  mart=ensembl)
b2_part2 <- getBM(attributes=c('ensembl_gene_id', 'description'),
                  values=a$X,
                  mart=ensembl)

anno_resdata_all<-merge(resdata_all,b2_part1,by.x="gene",by.y="ensembl_gene_id")

anno_resdata_p2<-merge(anno_resdata_df,b2_part2,by.x="gene",by.y="ensembl_gene_id")

library(dplyr)
anno_resdata <- distinct(anno_resdata)
View(anno_resdata)
consolidated_df <- anno_resdata %>%
  group_by(gene) %>%
  summarise(interpro_short_description = paste(unique(interpro_short_description), collapse = "; "))
View(consolidated_df)
anno_resdata_df<-merge(resdata,consolidated_df,by.x="gene",by.y="gene")
View(anno_resdata_df)

b2_part2 <- getBM(attributes=c('ensembl_gene_id', 'description'),
                  values=a$X,
                  mart=ensembl)

write.table(anno_resdata_p2,"deseq_control_v_Hh_iso",sep="\t")

res_2<-results(dds_deseq2, contrast=c("Treatment", "HhWT+isotype_control_mAb","HhWT+anti-p40"),alpha=0.05)

summary(res_2)

# out of 56576 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4430, 7.8%
# LFC < 0 (down)     : 3642, 6.4%
# outliers [1]       : 74, 0.13%
# low counts [2]     : 22568, 40%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
res_2= subset(res_2, padj<0.05)
res2 <- res_2[order(res_2$padj),]
res2_lfc5 <- subset(res2, log2FoldChange > 2 | log2FoldChange < -2)
View(res)
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds_deseq2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'gene'
View(resdata2)

anno_resdata2<-merge(resdata2,b2_part1,by.x="gene",by.y="ensembl_gene_id")

consolidated_df2 <- anno_resdata2 %>%
  group_by(gene) %>%
  summarise(interpro_short_description = paste(unique(interpro_short_description), collapse = "; "))
anno_resdata2_df <- distinct(anno_resdata2_df)
View(anno_resdata2_df)
anno_resdata2_df<-merge(anno_resdata2,consolidated_df2,by.x="gene",by.y="gene")
anno_resdata2_df <- anno_resdata2_df[, -c(26)]
anno_resdata2_df<-distinct(anno_resdata2_df)

write.table(anno_resdata2_df,"deseq_Hh_iso_v_Hh_p40_all",sep="\t")
