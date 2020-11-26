## Pre Decontamination Steps in QIIME2

# Separate NEC & AAT - get aat table
# Get taxonomy and phylo tree

#PREP
cd ~/Dropbox/AAT_microbiome/16S/
  source activate qiime2-2019.7
qiime --help
#separated metadata from seq lab into AAT_lab_metadata and NEC_lab_metadata
#other metadata file - AAT_workbook.txt

#Step 1 - Need to isolate AAT sequences from NEC
qiime feature-table filter-samples \
--i-table table.qza \
--m-metadata-file AAT_lab_metadata.txt \
--o-filtered-table aat-table.qza
qiime feature-table filter-seqs \
--i-data rep-seqs.qza \
--i-table aat-table.qza \
--o-filtered-data aat-rep-seqs.qza 

qiime feature-table summarize \
--i-table aat-table.qza \
--o-visualization aat-table.qzv \
--m-sample-metadata-file AAT_lab_metadata.txt
#237 samples, 14183 features, mean10736/sample, median 13342/sample
#from interactive table:
#highest control batch features 5210 - if cutoff at 5215 73%samples remain (~80% of AAT), 35% of features retained
#if lower end 500 - almost all AAT samples and ~control (93%total) samples remain, only 4% features retained
#midway: 2700, 81% of samples (almost all control, ~85%AAT) remain, 20% features retained
#2000 - 84%samples (90%AAT, 30% controls), 15% features
#1000 - 91%samples (almost all AAT, ~40%controls), 8%features
#10,000 - 63%samples, 60%features

qiime feature-table tabulate-seqs \
--i-data aat-rep-seqs.qza \
--o-visualization aat-rep-seqs.qzv
#Seq count - all 300/samples, total 14183

#Step 2 - Create a Phylogeny Tree
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences aat-rep-seqs.qza \
--o-alignment aligned-aat-rep-seqs.qza \
--o-masked-alignment masked-aligned-aat-rep-seqs.qza \
--o-tree unrooted-tree-aat.qza \
--o-rooted-tree rooted-tree-aat.qza

#Step 3 - Get Taxonomy
# Do after Training_ReferenceDatabase.R
qiime feature-classifier classify-sklearn \
--i-classifier gg-classifier.qza \
--i-reads aat-rep-seqs.qza \
--o-classification aat-gg-taxonomy.qza
qiime metadata tabulate \
--m-input-file aat-gg-taxonomy.qza \
--o-visualization aat-gg-taxonomy.qzv




