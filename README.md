# Host-microbiota-interaction

## RNA-seq basic steps
1. Trimmed with trimmomatic
2. Quality check with MultiQC
3. STAR two-step aligned : featureCounts output used
4. Differential expression determined with DESeq2

## 16 Analysis
V4 region of the microbiota processed with DADA2 method in QIIME2 using SILVA classifier

## Host-microbe interactions:
1. WGCNA
2. SparCC
3. Spearman correlation of WGCNA and SparCC
4. CCA
5. Random forest regression
