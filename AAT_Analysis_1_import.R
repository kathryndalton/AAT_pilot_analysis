#AAT Analysis: Part 1 Uploading Samples into QIIME and Demultiplexing 

#Prep
source activate qiime2-2019.7
qiime --help 
#Need Qiime1 to import dual_indexed files
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

cd ~/Downloads
md5 MacQIIME_1.9.1-20150604_OS10.7.tgz
tar -xvf MacQIIME_1.9.1-20150604_OS10.7.tgz
cd MacQIIME_1.9.1-20150604_OS10.7/
./install.s
source /macqiime/configs/bash_profile.txt
align_seqs.py -h ##to check if works

cd ~/Dropbox/AAT_microbiome/16S/
extract_barcodes.py --input_type barcode_paired_end -f forward.fastq -r reverse.fastq --bc1_len 12 --bc2_len 12
gzip barcodes.fastq

#Step 1 - Import and Demultiplex
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path for_upload \
  --output-path emp-paired-end-sequences.qza

qiime demux emp-paired \
--i-seqs emp-paired-end-sequences.qza \
--m-barcodes-file 20191023_Davis_run1_mapping_file.tsv \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--p-rev-comp-barcodes \
--o-per-sample-sequences demux-emp.qza \
--o-error-correction-details demux-details-emp.qza
##ALTERNATE###
qiime demux emp-paired \
--i-seqs emp-paired-end-sequences.qza \
--m-barcodes-file 20191023_Davis_run1_mapping_file.tsv \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--p-rev-comp-mapping-barcodes \
--o-per-sample-sequences demux-emp_alt.qza \
--o-error-correction-details demux-details-emp_alt.qza

qiime demux summarize \
--i-data demux-details-emp.qza \
--o-visualization demux-details.qzv
##288samples, range seq from 6-149512 med22115

qiime demux summarize \
--i-data demux-emp_alt.qza \
--o-visualization demux_alt.qzv
##288samples, range seq from 6-149512 med22115 - SAME, same samples TA1RV1 at bottom

qiime metadata tabulate \
--m-input-file demux-details-emp.qza \
--output-dir demux-details.qzv
## doesn't actually work
