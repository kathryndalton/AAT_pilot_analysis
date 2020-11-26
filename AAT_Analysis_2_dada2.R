#AAT Analysis: Part 2 DADA2/Quality Control & Denoise
# End result - feature table
#Contains both AAT and NEC sequences

#PREP
cd ~/Dropbox/AAT_microbiome/16S/
source activate qiime2-2019.7
qiime --help 

#Step 1 - Get metadata into QIIME
qiime metadata tabulate \
--m-input-file 20191023_Davis_run1_mapping_file.tsv \
--o-visualization tabulated-AAT_metadata.qzv

##Need to combine real metadata with this one from lab


#Step 2 - Denoise
# removes/corrects for "noisey" reads - removes low quality regions of the sequences
#based on results from Interactive Quality Plot from demux.qzc readout = 
#high quality scores for intial reads (keep all, trim 0), but start to decrease at end - kepp up to 300 (most)  
#Takes while to run ~10-20min
#ERROR RUNNING THIS IN R TERMINAL - NEED TO RUN IN SEPERATE TERMINAL WINDOW!! can go back to R for other commands
qiime dada2 denoise-single \
--i-demultiplexed-seqs demux-emp.qza \
--p-trim-left 0  \
--p-trunc-len 300  \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza  \
--o-denoising-stats stats-dada2.qza \
--verbose
#Removed TA001-RV1-Dust-Pre, TA001-RV1-Dost-Post, and NY91101-HV1-LR (both have very low sequences <20, <water)

#Step 3 - Tabulate stats table on filtering processes (which ones removed, etc)
qiime metadata tabulate \
--m-input-file stats-dada2.qza \
--o-visualization stats-dada2.qzv

#Step 4 - Generate Feature Tables
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file 20191023_Davis_run1_mapping_file.tsv
# total 285 recovered (3 filtered out, see abve) 11,178 mean seq per sample
#table.qza is everything - AAT and NEC

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
