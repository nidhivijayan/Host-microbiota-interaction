#!/bin/bash
#SBATCH --job-name=qiime2_1
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --account=xxx
#SBATCH --partition=standard
#SBATCH -o qiime2_1_%j.out
#SBATCH -e qiime2_1%j.err

conda init

conda activate qiime2-amplicon-2024.10

# qiime rescript get-silva-data \
#     --p-version '138.2' \
#     --p-target 'SSURef_NR99' \
#     --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
#     --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

# qiime rescript reverse-transcribe \
#     --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
#     --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

# #Removiing sequences based on size
# qiime rescript filter-seqs-length-by-taxon \
#     --i-sequences silva-138.2-ssu-nr99-seqs.qza \
#     --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
#     --p-labels Archaea Bacteria Eukaryota \
#     --p-min-lens 900 1200 1400 \
#     --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
#     --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza

# Extracting reference reads for 515f-806r (https://docs.qiime2.org/2024.10/tutorials/feature-classifier/)
#  qiime feature-classifier extract-reads \
#   --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza \
#   --p-f-primer GTGCCAGCMGCCGCGGTAA \
#   --p-r-primer GGACTACHVGGGTWTCTAAT \
#   --p-trunc-len 120 \
#   --p-min-length 100 \
#   --p-max-length 400 \
#   --o-reads silva-138.2-ref-seqs.qza

# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads /home/user/SILVA_db_qiime/silva-138-99-seqs-515-806.qza \
#   --i-reference-taxonomy /home/user/SILVA_db_qiime/silva-138-99-tax-515-806.qza \
#   --o-classifier silva-138-trained-classifier.qza

# For using clawabck method, I am extracting amplicon-specific reads, and then traiining the classifier with readytowear weighted reads 
# qiime feature-classifier extract-reads \
#     --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza \
#     --p-f-primer GTGYCAGCMGCCGCGGTAA \
#     --p-r-primer GGACTACNVGGGTWTCTAAT \
#     --p-n-jobs 2 \
#     --p-read-orientation 'forward' \
#     --o-reads silva-138.2-ssu-nr99-seqs-515f-806r.qza

# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads readytowear/data/gg_13_8/515f-806r/ref-seqs.qza \
#   --i-reference-taxonomy readytowear/data/gg_13_8/515f-806r/ref-tax.qza \
#   --i-class-weight readytowear/data/gg_13_8/515f-806r/soil-non-saline.qza \
#   --o-classifier gg138_v4_soil-non-saline_classifier.qza

# ALternatively, dereplicate with "uniq" mode, then extract amplicon reads, and then retrain the classifier. 

# echo -e 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath' > manifest2_STTR.tsv
# for FOR in /nfs/turbo/umms-youngvi/nidhi/STTR_16S/Young_M04695_R-443985636/part2/*_R1*fastq.gz;
# do
#  n=$(basename $FOR | cut -f1 -d_);
#  REV=${FOR/_R1_/_R2_}
#  echo $n
#  echo -e "$n\t$PWD/$FOR\t$PWD/$REV" >> /home/user/STTR_16S/manifest2_STTR.tsv;
# done

## Edit the file path in manifest.tsv to remove "/home/user"

# qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path /home/user/STTR_16S/manifest2_STTR.tsv \
#  --output-path sttr2-demux2.qza \
#  --input-format PairedEndFastqManifestPhred33V2

#   qiime demux summarize \
#  --i-data /home/user/STTR_16S/sttr2-demux2.qza \
#  --o-visualization /home/user/STTR_16S/demux2_2.qzv

#Removing the following samples because of low reads
#1511.cec	1972	1972
#1508.col	78	78

input=/home/user/STTR_16S/

# qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs $input/sttr2-demux2.qza \
#   --p-trunc-len-f 250 \
#  --p-trunc-len-r 230 \
#  --o-representative-sequences $input/sttr2-rep-seqs-dada2.qza \
#  --o-table $input/sttr2-table2-dada2.qza \
#  --o-denoising-stats $input/sttr2-stats2-dada2.qza \
#  --p-n-threads 5

# qiime feature-table merge \
# --i-tables sttr2-table2-dada2.qza \
# --i-tables table2-dada2.qza \
# --o-merged-table merged-tables.qza

# qiime feature-table merge-seqs \
# --i-data rep-seqs-dada2.qza \
# --i-data sttr2-rep-seqs-dada2.qza \
# --o-merged-data merged-rep-seqs.qza


 #picrust2_pipeline.py -s $input/exported-feature-seq/dna-sequences.fasta -i $input/exported-feature-table/feature-table.biom -o picrust2_out_pipeline -p 1

#add_description.py -i $input/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
#-o $input/picrust2_out_pipeline/EC_metagenome_out//pred_metagenome_unstrat_descrip.tsv.gz

# qiime metadata tabulate \
#   --m-input-file $input/merged-tables.qza \
#   --o-visualization $input/merged-tables.qzv

# qiime feature-table summarize \
#   --i-table $input/merged-tables.qza \
#   --o-visualization $input/merged-tables-summarized.qzv \
#   --m-sample-metadata-file $input/metadata-sttr.txt

# qiime feature-table tabulate-seqs \
#   --i-data merged-rep-seqs.qza \
#   --o-visualization merged-rep-seqs.qzv

# qiime phylogeny align-to-tree-mafft-fasttree \
#   --i-sequences merged-rep-seqs.qza \
#   --o-alignment aligned-sttr2-rep-seqs.qza \
#   --o-masked-alignment masked-aligned-sttr2-rep-seqs-dada2.qza \
#   --o-tree unrooted-tree-merged2.qza \
#   --o-rooted-tree rooted-tree-merged2.qza

# qiime diversity core-metrics-phylogenetic \
#   --i-phylogeny rooted-tree.qza \
#   --i-table filtered-table-d0.qza \
#   --p-sampling-depth 14090 \
#   --m-metadata-file 16s_metadata_d0.tsv \
#   --p-ignore-missing-samples \
#   --output-dir core-metrics-results-end

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-end/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file 16s_metadata_d0.tsv \
  --o-visualization core-metrics-results-end/weighted_unifrac_distance_matrix.qzv \
  --p-ignore-missing-samples \
  --p-pairwise \
  --m-metadata-column Treatment \

# qiime diversity beta-group-significance \
#   --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
#   --m-metadata-file metadata-sttr.txt \
#   --o-visualization core-metrics-results/weighted_UF_distance_matrix.qzv \
#   --p-pairwise \
#   --m-metadata-column treatment \

#   qiime diversity beta-group-significance \
#   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
#   --m-metadata-file metadata-sttr.txt\
#   --o-visualization core-metrics-results/unweighted_UF_distance_matrix.qzv \
#   --p-pairwise \
#   --m-metadata-column treatment \

# qiime tools export \
# --input-path rooted_tree.qza \
# # --output-path tree

# classifier=/home/user/SILVA_db_qiime/silva-138-trained-classifier.qza
# qiime feature-classifier classify-sklearn \
#   --i-classifier $classifier \
#   --i-reads merged-rep-seqs.qza \
#   --o-classification merged-taxonomy-silva138-1.qza

# qiime taxa barplot \
#   --i-table table2-dada2.qza \
#   --i-taxonomy taxonomy-silva138-1.qza\
#   --m-metadata-file metadata-sttr.txt \
#   --o-visualization taxa-bar-plots-silva138-ed.qzv

# qiime feature-classifier classify-sklearn \
#   --i-classifier /home/user/qiime2/silva-138-trained-classifier.qza \
#   --i-reads /home/user/STTR_16S/exported_rep/rep-seqs-ed.qza \
#   --o-classification silva1-taxonomy-ed.qza


# qiime feature-table filter-samples \
#     --i-table table-dada2.qza \
#     --m-metadata-file 16s_metadata_d0.tsv \
#     --o-filtered-table filtered-table-d0.qza

