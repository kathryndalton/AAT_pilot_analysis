AAT Analysis Part 1 - QIIME Processing
================

Bioinformatic Processing Steps in QIIME
=======================================

Need to start with fastq files (etiehr stored on local computer or downloaded from NCBI SRA) \# Notes these files contain both AAT and NEC files Run all qiime files in command line (Terminal, etc) All .qza and .qzv files can be viewed on view.qiime2.org (.qza artifacts, .qzv visualization)

It is helpful to review background on QIIME before starting to understand the overall proces flow, see <https://docs.qiime2.org/2020.8/tutorials/moving-pictures/>

Overview of Steps
-----------------

#### 1. Import files into QIIME

#### 2. Demultiplex using DADA2 (match sequences to samples)

#### 3. Denoise using DADA2 (remove low quality reads)

#### 4. Seperate AAT samples from NEC samples (or other study samples)

#### 5. Match to Phylogeny

#### 6. Match to Taxonomy

#### Optional Step 7: Manually BLAST Staphlococcus samples to species level

End Products
------------

#### Feature Table

#### Taxonomy Table

#### Phylogeny Table

Step 1. Import files into QIIME
-------------------------------

Need to download QIIME (both 1 and 2) before proceeding. See <https://docs.qiime2.org/2020.8/install/> <http://qiime.org/install/install.html>

Will need to change "Downloads" and "Dropbox" to whereever you saved your files Also make sure names of files (barcode\_paired\_end, etc) matches files name

QIIME 1 is needed since are sequencing files are dual\_indexed, meaning the forward and reverse sequencing primers have different barcodes and we need to uniquely match based on BOTH of these barcode. QIIME2 is apparently working on incorporating this feature - will update this file if they do.

``` bash
############################################
#Prep
source activate qiime2-2019.7
qiime --help 
#Need Qiime1 to import dual_indexed files
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
cd ~/Downloads   ### change to wherever QIIME1 files are located
md5 MacQIIME_1.9.1-20150604_OS10.7.tgz
tar -xvf MacQIIME_1.9.1-20150604_OS10.7.tgz
cd MacQIIME_1.9.1-20150604_OS10.7/
./install.s
source /macqiime/configs/bash_profile.txt
align_seqs.py -h ##to check if works

cd ~/Dropbox/AAT_microbiome/16S/ ## change to wherever your sequencing files (and barcodes) are located
extract_barcodes.py --input_type barcode_paired_end -f forward.fastq -r reverse.fastq --bc1_len 12 --bc2_len 12
gzip barcodes.fastq
#################################################

#Step 1 - Import sequences
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path for_upload \
  --output-path emp-paired-end-sequences.qza
  
# Import metatable file (for now start with one just from lab - will combine full/extended one later in R)
  qiime metadata tabulate \
--m-input-file 20191023_Davis_run1_mapping_file.tsv \  ## rename file to match - make sure this is in working directory
--o-visualization tabulated-AAT_metadata.qzv
```

Step 2. Demultiplexing into QIIME
---------------------------------

Demultiplexing is the process of matching the sequences to our samples based on barcodes added during the sequencing process. For this, we need not only the sequences, but the metadata from the sequencing facility which will have the unique barcodes for every sample. Our samples are paired-end reads, meaning there are forward and reads per sample.

For quality control, we'll generate a summary of the demultiplexing results. This allows you to determine how many sequences were obtained per sample, and also to get a summary of the distribution of sequence qualities at each position in your sequence data.

``` bash

## Step 2: Demultiplexing

qiime demux emp-paired \
--i-seqs emp-paired-end-sequences.qza \
--m-barcodes-file 20191023_Davis_run1_mapping_file.tsv \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--p-rev-comp-barcodes \
--o-per-sample-sequences demux-emp.qza \
--o-error-correction-details demux-details-emp.qza

##ALTERNATE###  --- gives same results, recommend to use both ways with your data to confirm the same (use demux summarize to confirm)
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
```

Step 3. Denoise using DADA2
---------------------------

Removes/corrects for "noisey" reads - removes low quality regions of the sequences, chimeras, etcs THis is based on results from Interactive Quality Plot from demux.qzc readout - need to explore and create your own cut-off point --p-trim-left m, which trims off the first m bases of each sequence (based on quality of initial few reads) --p-trunc-len n which truncates each sequence at position n, based on where quality starts to drop off For this dataset, there was high quality scores for intial reads (keep all, trim 0), but start to decrease at end - kepp up to 300 (most)

#### Takes while to run ~10-20min

ERROR RUNNING THIS IN R TERMINAL - NEED TO RUN IN SEPERATE TERMINAL WINDOW!! can go back to R for other commands

#### End result = Feature Table

``` bash

# Step 3: Denoise
qiime dada2 denoise-single \
--i-demultiplexed-seqs demux-emp.qza \
--p-trim-left 0  \
--p-trunc-len 300  \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza  \
--o-denoising-stats stats-dada2.qza \
--verbose
#Removed TA001-RV1-Dust-Pre, TA001-RV1-Dost-Post, and NY91101-HV1-LR (both have very low sequences <20, <water)

#Tabulate stats table on filtering processes (which ones removed, etc)
qiime metadata tabulate \
--m-input-file stats-dada2.qza \
--o-visualization stats-dada2.qzv

#Evaulate Feature Tables and Final Sequences
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file 20191023_Davis_run1_mapping_file.tsv
# total 285 recovered (3 filtered out, see abve) 11,178 mean seq per sample
#table.qza is everything - AAT and NEC
qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
```

Step 4. Seperate AAT samples from NEC samples (or other study samples)
----------------------------------------------------------------------

Currently these sequences and feature table combine both samples from AAT with NEC (run togehter on sample batch / plates). Need to seperate them out before prcoeeding further with denoising and analysis. We are filtering samples by the full metatable table, where any samples not listed in this file are removed.

``` bash

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
```

Step 5. Match to Phylogeny
--------------------------

Needed for any phylogenic -based metric analysis (Faith's, weight UniFrac, etc). The pipeline uses the mafft program to perform a multiple sequence alignment of the sequences in our FeatureData\[Sequence\] to create a FeatureData\[AlignedSequence\] QIIME 2 artifact. Next, the pipeline masks (or filters) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree. Following that, the pipeline applies FastTree to generate a phylogenetic tree from the masked alignment. The FastTree program creates an unrooted tree, so in the final step in this section midpoint rooting is applied to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree.

``` bash
 
 # Phylogeny
 qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences aat-rep-seqs.qza \
--o-alignment aligned-aat-rep-seqs.qza \
--o-masked-alignment masked-aligned-aat-rep-seqs.qza \
--o-tree unrooted-tree-aat.qza \
--o-rooted-tree rooted-tree-aat.qza
```

Step 6. Match to Taxonomy
-------------------------

We need to assign taxonomy to the sequences in our feature table. Can use an already trained (created) classifer (reference databse), but it is better to "train" one on your dataset. \#\#\#\# This is a multi-step process which take over a day to run!

Before you run this, you need to download the full taxonomy reference database either from Greengenes or SILVA (or other site) - I've found Greengenes better at detecting Staph species. Do this by going to <https://docs.qiime2.org/2020.8/data-resources/> - Go to Marker Gene Reference Databases -&gt; Greenegenes (16S rRNA) - click on 13\_8 (or whatever id the most recent) \*\* apparently was issues with 13\_8, so GG went back to 13\_5) - unzip and will have 5 folders (otus, rep\_set, rep\_set\_aligned, trees, and taxonomy, plus a notes file) OR <https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/> - click on most recent (one with highest number) - download gg\_\#\#\#\_otus.tar.gz (whatever number) - will have same 5 folders

We'll be using 99\_ files = reference sequences clustered at 99% sequence similarity

##### Can also use SILVA as comparison - google database

The we filter the reference database according to our sequencing protocols - this case 300bp V1-3 using 27F/534R primer pair 27F = AGAGTTTGATCMTGGCTCAG 534R = ATTACCGCGGCTGCTGG May want to confirm with sequencing facility what primers they used

``` bash
# Creating a Pre-Trained Classifier ##
## import downloaded reference files
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path 99_otus.fasta \
--output-path 99_otus.qza
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path 99_otu_taxonomy.txt \
--output-path ref-taxonomy.qza
# Extract Reference Reads
# sequence reads that we’re trying to classify are 300-base single-end reads that were amplified with the 27F/534R primer pair for 16S rRNA gene sequences V1-3
# optimize for that here by extracting reads from the reference database based on matches to this primer pair, and then slicing the result to 300 bases (-p-trunc-len)
# min & max length exclude simulated amplicons that are far outside of the anticipated length distribution using those primers / non-targets
qiime feature-classifier extract-reads \
--i-sequences 99_otus.qza \
--p-f-primer AGAGTTTGATCCTGGCTCAG \
--p-r-primer ATTACCGCGGCTGCTGG \
--p-trunc-len 300 \
--p-min-length 200 \
--p-max-length 500 \
--o-reads ref-seqs.qza

## Train Classifier 
## THIS TAKES A LONG TIME ~ ~ About one day or overnight depending on your computer
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier gg-classifier.qza

#Test classifier with your own data - Assign Taxonomy 
#In this case will use AAT rep sequences
qiime feature-classifier classify-sklearn \
--i-classifier gg-classifier.qza \
--i-reads aat-rep-seqs.qza \
--o-classification aat-gg-taxonomy.qza
qiime metadata tabulate \
--m-input-file aat-gg-taxonomy.qza \
--o-visualization aat-gg-taxonomy.qzv
```

Now you have all your three main files needed for analysis
----------------------------------------------------------

#### 1. Feature Table (list of taxa numbers by sample)

#### 2. Taxonomy Table (what bacteria matches what taxa)

#### 3. Phylogeny Table (how taxa are related to eachother phylogenetically)

Next STEP is moving these files into R for further processing and analysis
--------------------------------------------------------------------------
