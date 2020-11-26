# Decontaminating - removing likely contaminants
#used to ID contaminant sequences based on frequency of each sequence/OTU/ASV 
#in the input feature tables as a function of the concentration of amplified DNA
#Needs: feature table, concentration DNA, negative controls
#Method - based on either conc or negs, or combined(using Fishers)
#ID contaminants from sample data (con/negs) then go back to feature table to remove
#Best way to do it is to start with a phyloseq object - a file that combines otu feature table + sample data + tax data

#Download & Set Up#
#Phyloseq
setwd("~/Dropbox/AAT_microbiome/16S")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("phyloseq")
packageVersion('phyloseq')
library(phyloseq)
#Decontam
BiocManager::install("decontam")
packageVersion("decontam")
library(decontam)
seqtab<-read_qza("aat-table.qza")
tree<-read_qza("rooted-tree-aat.qza")
taxonomy<-read_qza("aat-gg-taxonomy.qza")
taxtable<-taxonomy$data %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
aat_lab_metadata <-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/AAT_lab_metadata.txt") 
aat_lab_metadata$X16S_qPCR_copy_num_per_ul_DNA<-as.numeric(aat_lab_metadata$X16S_qPCR_copy_num_per_ul_DNA)

#Create Phyloseq Object = combination all files (otu table, taxonomy, metadata)
##DO NOT RUN 
#phylo<-qza_to_phyloseq("aat-table.qza", "rooted-tree-aat.qza", "aat-taxonomy.qza","AAT_lab_metadata.txt")
#This version doesn't work - DROPS OFF FIRST SAMPLE (control.batch10) -> SEE BELOW FOR TROUBLESHOOTING
phylo1<-phyloseq(
  otu_table(seqtab$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
  sample_data(aat_lab_metadata %>% as.data.frame() %>% column_to_rownames("#SampleID"))
)
phylo1
dim(seqtab$data)
dim(aat_lab_metadata)
sample_names(phylo1)
table(sample_data(phylo1)$batch)
#table missing from DADA2 - #Removed TA001-RV1-Dust-Pre, TA001-RV1-Dost-Post

#Simple Phyloseq Analysis
plot_bar(phylo1, fill = "Family")
plot_tree(phylo1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_heatmap(phylo1)
plot_heatmap(phylo1, taxa.label="Phylum")

#Inspect Library Size
df <- as.data.frame(sample_data(phylo1)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phylo1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Study)) + 
  geom_point()
ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType)) + 
  geom_point()
#Index = ordered based on library size


#####################################################################################

#### FINAL STEPS FOR HOW TO HANDLE CONTAMINANTS

##Step 1 - ID contam using combo batch for lab controls + DNA conc data               166 contam
##Step 2 - ID contam using combo batch for lab controls + qPCR conc data              149 contam
##Step 3 - Combine steps 2 + 3 - ID contaminants found on EITHER methods              261 contam (54 shaared)
##Step 4 - Remove those contaminants from phyloseq object                             13922 taxa remain
##Step 5 - ID contaminants in new phyloseq objects using field blanks - loose ID      188 contam
##Step 6 - Remove field contaminants to get final non-contaminated phyloseq object    13734 taxa remian 
# phylo.final  = final object
## Make list of contaminants at each step!!
## NEED TO FIGURE OUT HOW TO SAVE phylo.final AS ACTUAL FILE SO DON'T HAVE TO RERUN WHEN CLEARING ENV 

##Step 1 - ID contam using combo batch for lab controls + DNA conc data 
sample_data(phylo1)$is.neg <- sample_data(phylo1)$Study=="control"
contamdf.combobatch <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", 
                                     batch = "batch")
table(contamdf.combobatch$contaminant)  #using dna 166
head(which(contamdf.combobatch$contaminant)) # 832nd most abundant
taxtable1<-taxtable %>% select(-Confidence) %>% column_to_rownames("Feature.ID")
contam_asvs1 <- row.names(contamdf.combobatch[contamdf.combobatch$contaminant == TRUE, ])
contam_taxa1<-taxtable1[row.names(taxtable1) %in% contam_asvs1, ]
sort(table(contam_taxa1$Genus), decreasing = TRUE)

##Step 2 - ID contam using combo batch for lab controls + qPCR conc data
phylo_qpcr <- prune_samples(!is.na(sample_data(phylo1)$X16S_qPCR_copy_num_per_ul_DNA), phylo1)
phylo_qpcr #should have 233 samples (-3 NA samples)
phylo1
contamdf.combobatch.qpcr <- isContaminant(phylo_qpcr, method="combined", neg="is.neg", conc="X16S_qPCR_copy_num_per_ul_DNA", 
                                          batch = "batch")
table(contamdf.combobatch.qpcr$contaminant) #using qpcr 149
head(which(contamdf.combobatch.qpcr$contaminant))  #591st most abundant
contam_asvs2 <- row.names(contamdf.combobatch.qpcr[contamdf.combobatch.qpcr$contaminant == TRUE, ])
contam_taxa2<-taxtable1[row.names(taxtable1) %in% contam_asvs2, ]
sort(table(contam_taxa2$Genus), decreasing = TRUE)

##Step 3 - Combine steps 2 + 3 - ID contaminants found on EITHER methods 
c1<-rownames_to_column(contamdf.combobatch, var="rowname")
c2<-rownames_to_column(contamdf.combobatch.qpcr, var="rowname")
c1$contaminant.d<-c1$contaminant
c2$contaminant.q<-c2$contaminant
c3<-full_join(c1,c2, by="rowname") #x=dna, y=qpcr  eg. c3$contaminant.x = c3$contaminant.d
table(c3$contaminant.d)
table(c3$contaminant.q) 
c3$contaminant.both<-ifelse(c3$contaminant.d==TRUE&c3$contaminant.q==TRUE, TRUE, FALSE)  
c3$contaminant.either<-ifelse(c3$contaminant.d==TRUE, TRUE, 
                              ifelse(c3$contaminant.q==TRUE, TRUE, FALSE))  
table(c3$contaminant.both)  # 54 shared between DNA + qPCR results
table(c3$contaminant.either)  #261 total contaminants IDed using either conc method

##Step 4 - Remove those contaminants from phyloseq object 
phylo1
ps.noncontam.labcontconc <- prune_taxa(!c3$contaminant.either, phylo1)
ps.noncontam.labcontconc # taxa should be less (original 14183 - 261above = 13922 total taxa)

##Step 5 - ID contaminants in new phyloseq objects using field blanks
sample_data(ps.noncontam.labcontconc)$field.control <- sample_data(ps.noncontam.labcontconc)$SampleType=="Environmental Blank"
table(sample_data(ps.noncontam.labcontconc)$field.control)  #Missing TA001 RV1 Blank == 12 total blanks
contamdf.field <- isContaminant(ps.noncontam.labcontconc, method="prevalence", neg="field.control")
table(contamdf.field$contaminant) #188 contaminates IDed
head(which(contamdf.field$contaminant)) #705th most abundant is top one
contam_asvs3 <- row.names(contamdf.field[contamdf.field$contaminant == TRUE, ])
contam_taxa3<-taxtable1[row.names(taxtable1) %in% contam_asvs3, ]
sort(table(contam_taxa2$Genus), decreasing = TRUE)

##Step 6 - Remove field contaminants to get final non-contaminated phyloseq object
phylo1   #started with 14183 taxa
ps.noncontam.labcontconc  #went down to 13922
phylo.final <- prune_taxa(!contamdf.field$contaminant, ps.noncontam.labcontconc)
phylo.final # final tally 13734   (13922 - 188)


##################################################################################################################


###  Background Playing with Decontam and Checking Data ####


#Identify Contaminants based on Frequency (DNA Concentration)
# X16S_final_library_concentration_ng_ul = sample variable that holds the concentration information:
# = quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
##All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
contamdf.freq <- isContaminant(phylo1, method="frequency", conc="X16S_final_library_concentration_ng_ul")
head(contamdf.freq)
# output = 
# $p which containts the probability that was used for classifying contaminants
# $contaminant which contains TRUE/FALSE classification values with TRUE indicating that the statistical 
# evidence that the associated sequence feature is a contaminant exceeds the user-settable threshold. 
# As we did not specify the threshold, the default value of threshold = 0.1 was used, 
#and $contaminant=TRUE if $p < 0.1
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))
# 320 out of 14183 (~2%) of ASVs were classified as contaminants based on DNA concentration
# This includes the 41st & 56st most abundant ASVs
plot_frequency(phylo1, taxa_names(phylo1)[c(2,41)], conc="X16S_final_library_concentration_ng_ul") + 
  xlab("X16S final library concentration ng/ul")
# shows frequency is expected to be inversely proportional to input DNA concentration, as 
# contaminating DNA will make up a larger fraction of the total DNA in samples with very little total DNA
# graph on right is supposed to fit model (contaminant), while one of left does not
set.seed(100)
plot_frequency(phylo1, taxa_names(phylo1)[sample(which(contamdf.freq$contaminant),3)], conc="X16S_final_library_concentration_ng_ul") +
  xlab("X16S final library concentration ng/ul")
#graphs of other contaminants
# Remove ID contaminants from the phyloseq object
phylo1
ps.noncontam.freq <- prune_taxa(!contamdf.freq$contaminant, phylo1)
ps.noncontam.freq # taxa should be less and should match FALSE for table(contamdf.freq$contaminant) above

#Identify Contaminants based on Prevalence - Negative Controls
# Study is the sample variable that holds the negative control information, need to convert to logical
sample_data(phylo1)$is.neg <- sample_data(phylo1)$Study=="control"
phylo1 #should now have 47 sample variables (+1)
contamdf.prev <- isContaminant(phylo1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
# prev has less number of contaminants (41 vs 320 for dna) and less abdundant (408 highest)
#the default threshold for a contaminant is that it reaches a probability of 0.1 in Fishers test
#revalence test there is a special value worth knowing, threshold=0.5, that will identify as contaminants 
#all sequences thare are more prevalent in negative controls than in positive samples
contamdf.prev05 <- isContaminant(phylo1, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))
ps.pa <- transform_sample_counts(phylo1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Study == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Study == "AAT", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#shows number of times several of these taxa were observed in negative controls and positive samples
phylo1
ps.noncontam.prev <- prune_taxa(!contamdf.prev$contaminant, phylo1)
ps.noncontam.prev # taxa should be less and should match FALSE for table(contamdf.freq$contaminant) above


#Other Methods to ID Contaminants
# combined = freq & prev probabilities are combined by Fishers and used to ID
contamdf.combo <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul")
table(contamdf.combo$contaminant)
head(which(contamdf.combo$contaminant))
contamdf.combo05 <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", threshold=0.5)
table(contamdf.combo05$contaminant)
head(which(contamdf.combo05$contaminant))
# minimum = min of the freq + prev is used to ID 
contamdf.min <- isContaminant(phylo1, method="minimum", neg="is.neg", conc="X16S_final_library_concentration_ng_ul")
table(contamdf.min$contaminant)
head(which(contamdf.min$contaminant))
# either = contaminants are called if IDd by EITHER freq or prev 
contamdf.either <- isContaminant(phylo1, method="either", neg="is.neg", conc="X16S_final_library_concentration_ng_ul")
table(contamdf.either$contaminant)
head(which(contamdf.either$contaminant))
# both = contaminants are called if IDd by BOTH freq and prev
contamdf.both <- isContaminant(phylo1, method="both", neg="is.neg", conc="X16S_final_library_concentration_ng_ul")
table(contamdf.both$contaminant)
head(which(contamdf.both$contaminant))
#End Results - # of contaminants removed / Highest Abundant
# Freq/DNA conc = 320  /  41
# Prev = 41   /  408
# Prev@05 = 126  /  200
# Combo = 122  / 56   **** think this is the best option moving forward
# Combo@05  = 872 / 10
# Min = 309  / 41
# Either = 343  / 41
# Both = 2  / 1549
#Checking combination results
ps.pa1 <- transform_sample_counts(phylo1, function(abund) 1*(abund>0))
ps.pa.neg1 <- prune_samples(sample_data(ps.pa1)$Study == "control", ps.pa1)
ps.pa.pos1 <- prune_samples(sample_data(ps.pa1)$Study == "AAT", ps.pa1)
df.pa1 <- data.frame(pa.pos1=taxa_sums(ps.pa.pos1), pa.neg1=taxa_sums(ps.pa.neg1),
                    contaminant=contamdf.combo$contaminant)
ggplot(data=df.pa1, aes(x=pa.neg1, y=pa.pos1, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Other options 
#batch = seq run, test done independently, can be seq run, plate or extract batch
contamdf.combobatch <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", 
                                batch = "batch")
table(contamdf.combobatch$contaminant)
head(which(contamdf.combobatch$contaminant))
contamdf.combobatch05 <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", 
                                     threshold = 0.5, batch = "batch")
table(contamdf.combobatch05$contaminant)
head(which(contamdf.combobatch05$contaminant))
# Results  166 removed / 832 highest abundant, for threshold 0.5 = 1060 / 13
#batch.combine = prob each batch combined to create new threshold, one above is minimum
contamdf.combobatchprod <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", 
                                        batch = "batch", batch.combine = "product")
table(contamdf.combobatchprod$contaminant)
head(which(contamdf.combobatchprod$contaminant))  # 504 removed, 21st abundant
contamdf.combobatchfish <- isContaminant(phylo1, method="combined", neg="is.neg", conc="X16S_final_library_concentration_ng_ul", 
                                        batch = "batch", batch.combine = "fisher")
table(contamdf.combobatchfish$contaminant)
head(which(contamdf.combobatchfish$contaminant)) # 11 removed, 1080st abundant 
# Final Check - using combo + batch  (default minimum)
ps.pa2 <- transform_sample_counts(phylo1, function(abund) 1*(abund>0))
ps.pa.neg2 <- prune_samples(sample_data(ps.pa2)$Study == "control", ps.pa2)
ps.pa.pos2 <- prune_samples(sample_data(ps.pa2)$Study == "AAT", ps.pa2)
df.pa2 <- data.frame(pa.pos2=taxa_sums(ps.pa.pos2), pa.neg2=taxa_sums(ps.pa.neg2),
                    contaminant=contamdf.combobatch$contaminant)
ggplot(data=df.pa2, aes(x=pa.neg2, y=pa.pos2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#     Do both qpcr (from raw samples) and post-pcr amplification dna conc 
phylo_qpcr <- prune_samples(!is.na(sample_data(phylo1)$X16S_qPCR_copy_num_per_ul_DNA), phylo1)
phylo_qpcr #should have 233 samples (-3 NA samples)
phylo1
contamdf.freq.qpcr<- isContaminant(phylo_qpcr, method="frequency", conc="X16S_qPCR_copy_num_per_ul_DNA")
table(contamdf.freq.qpcr$contaminant) #using qpcr 166 removed
table(contamdf.freq$contaminant) #using dna conc = 320 removed

plot_frequency(phylo1, taxa_names(phylo1)[sample(which(contamdf.freq$contaminant),3)], conc="X16S_final_library_concentration_ng_ul") +
  xlab("X16S final library concentration ng/ul")
plot_frequency(phylo_qpcr, taxa_names(phylo_qpcr)[sample(which(contamdf.freq.qpcr$contaminant),3)], conc="X16S_qPCR_copy_num_per_ul_DNA") +
  xlab("X16S qPCR copies/ul DNA")

contamdf.combobatch.qpcr <- isContaminant(phylo_qpcr, method="combined", neg="is.neg", conc="X16S_qPCR_copy_num_per_ul_DNA", 
                                     batch = "batch")
table(contamdf.combobatch.qpcr$contaminant) #using qpcr 149
table(contamdf.combobatch$contaminant) #using dna 166
#      Keep contaminant if found on 1) both runs or 2) either run **see total numbers to decide** 
c1<-rownames_to_column(contamdf.combobatch, var="rowname")
c2<-rownames_to_column(contamdf.combobatch.qpcr, var="rowname")
c1$contaminant.d<-c1$contaminant
c2$contaminant.q<-c2$contaminant
c3<-full_join(c1,c2, by="rowname") #x=dna, y=qpcr  eg. c3$contaminant.x = c3$contaminant.d
table(c3$contaminant.d)
table(c3$contaminant.q) 
c3$contaminant.both<-ifelse(c3$contaminant.d==TRUE&c3$contaminant.q==TRUE, TRUE, FALSE)  
c3$contaminant.either<-ifelse(c3$contaminant.d==TRUE, TRUE, 
                              ifelse(c3$contaminant.q==TRUE, TRUE, FALSE))  
table(c3$contaminant.both)
table(c3$contaminant.either)


ps.pa <- transform_sample_counts(phylo1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Study == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Study == "AAT", ps.pa)
df.pa3 <- data.frame(pa.pos2=taxa_sums(ps.pa.pos2), pa.neg2=taxa_sums(ps.pa.neg2),
                     contaminant=c3$contaminant.either)
ggplot(data=df.pa2, aes(x=pa.neg2, y=pa.pos2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
hist(c3$p.x) #Confirms 0.5 threshold
hist(c3$p.y)

# Remove ID contaminants from the phyloseq object
phylo1
ps.noncontam.labcontconc <- prune_taxa(!c3$contaminant.either, phylo1)
ps.noncontam.labcontconc # taxa should be less (original 14183 - 261above = 13922 total taxa)




#Also figure out how to differentiate extract vs seq vs field neg controls
#RESULT LITTLE DIFFERENCE BETWEEN SEQ AND EXTRACT CONTROLS SEPERATELY VS LAB CONTROLS COMBINED
#Extract controls is primary driver of lab controls - seq controls has very little contamintion 

#METHODS TO DO - NOT NEEDED FOR ACTUAL FINAL ANALYSIS
#     rerun is.cont for each type of control, using different parameters for each, less strict field
#        look to see what and how much you are losing at each step 
sample_data(phylo1)$seq.control <- sample_data(phylo1)$batch=="seq"
sample_data(phylo1)$extract.control <- sample_data(phylo1)$SampleType=="Extract Control"
sample_data(phylo1)$field.control <- sample_data(phylo1)$SampleType=="Environmental Blank"
phylo1 #should now have 50 sample variables (47+3)
#Stage 1 - Strict using seq controls
#Stage 2 - Strict using extract controls
#Stage 3 - Conservative using field controls

#play around with when to remove based on dna/qpcr

contamdf.prev11 <- isContaminant(phylo1, method="prevalence", neg="seq.control")
contamdf.prev15 <- isContaminant(phylo1, method="prevalence", neg="seq.control", threshold=0.5)
table(contamdf.prev11$contaminant) # 4 contaminants
table(contamdf.prev15$contaminant)  #12 contaminants 
contamdf.prev21 <- isContaminant(phylo1, method="prevalence", neg="extract.control")
contamdf.prev25 <- isContaminant(phylo1, method="prevalence", neg="extract.control", threshold=0.5)
table(contamdf.prev21$contaminant) # 78 contaminants
table(contamdf.prev25$contaminant)  #142 contaminants 
table(contamdf.prev$contaminant) #using both as "lab controls" 41
table(contamdf.prev05$contaminant)  # 126 
contam_ext <- row.names(contamdf.prev25[contamdf.prev25$contaminant == TRUE, ])
contam_lab <- row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
shared_control<-contam_ext[contam_ext %in% contam_lab] #106 shared contaminants bet using just extract controls vs both lab controls 
cont1<-c(contam_ext, contam_lab)
cont_either<-unique(cont1) #162 total (142+126 - 106)   so 20 unique from lab, 36 from extract  
contam_seq <- row.names(contamdf.prev15[contamdf.prev15$contaminant == TRUE, ])
shared_control1<-contam_seq[contam_seq %in% contam_ext] # 3 contaminants shared between seq controls $ extract controls
shared_control2<-contam_seq[contam_seq %in% contam_lab] #8 shared bet seq controls & lab controls  5 are unique just to lab&seq, 3 are shared all 3
shared_control3<-contam_seq[contam_seq %in% shared_control]
cont2<-c(contam_ext, contam_lab, contam_seq)
cont_any<-unique(cont2) #166 total (142+126+12 - 106(lab-ext) - 3(all 3) -5(lab-seq) ) 
s1<-unique(c(shared_control1, shared_control2, shared_control3)) #shows 8 total found within lab
cont3<-c(contam_ext, contam_seq)
cont_s.e<-unique(cont3) #151 unique contam betw seq & extract controls >> 126 from lab controls combined (139 ext + 9 seq + 3 shared OR 142+12-3)
###Overall recomm same strict policy for both seq&ext controls, but doing is.contaminants sequenciallly not using combined lab controls
s1<-unique(c(shared_control1, shared_control2, shared_control3))
df.pa3 <- data.frame(pa.pos2=taxa_sums(ps.pa.pos2), pa.neg2=taxa_sums(ps.pa.neg2),
                     contaminant=contamdf.combobatch$contaminant)
ggplot(data=df.pa2, aes(x=pa.neg2, y=pa.pos2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Need to keep data.set with ASV as rownames (use contamf.prev15 and 25)
#What I need to do is figure out how to just have two contaminant rows - one for contam.seq T/F, one for contam.ext T/F 
#then create seonc for contam.both or contam.either
q1<-rownames_to_column(contamdf.prev15, var="rowname")
q2<-rownames_to_column(contamdf.prev25, var="rowname")
q1$contaminant.s<-q1$contaminant
q2$contaminant.e<-q2$contaminant
q4<-full_join(q1,q2, by="rowname") #x=seq, y=ext  eg. q5$contaminant.s = q5$contaminant.x
table(q4$contaminant.s)
table(q4$contaminant.e) 
q4$contaminant.both<-ifelse(q4$contaminant.s==TRUE&q4$contaminant.e==TRUE, TRUE, FALSE)  
q4$contaminant.either<-ifelse(q4$contaminant.s==TRUE, TRUE, 
                              ifelse(q4$contaminant.e==TRUE, TRUE, FALSE))  
table(q4$contaminant.both)
table(q4$contaminant.either)
df.pa3 <- data.frame(pa.pos2=taxa_sums(ps.pa.pos2), pa.neg2=taxa_sums(ps.pa.neg2),
                     contaminant=q4$contaminant.either)
ggplot(data=df.pa2, aes(x=pa.neg2, y=pa.pos2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
hist(q4$p.x) #Confirms 0.5 threshold
hist(q4$p.y)
# Remove ID contaminants from the phyloseq object
phylo1
ps.noncontam.freq <- prune_taxa(!contamdf.freq$contaminant, phylo1)
ps.noncontam.freq # taxa should be less and should match FALSE for table(contamdf.freq$contaminant) above

#Next Steps: When to remove based on DNA/qPCR contam, when to remove based on field blanks 
#Step 1 - IDed contamin based on seq and extract blanks SOLUTION - do seperately than EITHER combine 151 total
#### NEED TO REDO USING BATCHES!!
#Step 2 - combine these contaminants with comtaminants found using DNA + qPCR either
cb1<-c(contam_cb, contam_cb_qpcr)
cb_either<-unique(cb1) #261 total (166+149 - 54)               
either_tax_cb<-taxtable1[row.names(taxtable1) %in% cb_either, ] 








#Checking appropriate p value selection 
# simplest and most useful evaluation of method is to inspect the distribution of scores assigned by the method.
#Expectation is that there will be a strong mode at low scores. In the cleanest cases the 
#distribution will be clearly bimodal, while in other datasets the high-score mode is more wide and 
#diffuse. However, the low-score mode should be there, and should be used to set the P* score threshold 
#for identifying contaminants
hist(contamdf.combobatch$p, 100) #slight bimodal
hist(contamdf.combo$p, 100) #more bimodal
hist(contamdf.prev$p, 100) # these two look much better
hist(contamdf.prev05$p, 100)
hist(contamdf.combobatch05$p, 100) #this one looks worse

#Added Steps - get list of contaminant
tax_table(phylo1)[contamdf.combobatch$contaminant,] # to get genera corresponding to contaminants
ps.contam <- prune_taxa(contamdf.combobatch$contaminant, phylo1) # create a phyloseq objects w just to contaminants
head(tax_table(ps.contam))
ps.contam_list<-as.data.frame(tax_table(ps.contam)) 
#Way 2
contam <- row.names(contamdf.combobatch[contamdf.combobatch$contaminant == TRUE, ])
taxtable1<-taxtable %>% column_to_rownames("Feature.ID")
contam_list<-taxtable1[row.names(taxtable1) %in% contam, ]  %>% select(-Confidence)
head(contam_list)





#TROUBLESHOOTING QZA_TO_PHYLOSEQ
#table missing from DADA2 - #Removed TA001-RV1-Dust-Pre, TA001-RV1-Dost-Post
sample_data<- as.data.frame(sample_data(phylo))
table(sample_data$control)
table(aat_lab_metadata$Study)
dim(seqtab$data)
dim(aat_lab_metadata)
head(sample_data(phylo))
head(otu_table(phylo))
sample_names(phylo)
#for some reason control.batch.10 is cut off on merging into phyloseq object
# option is individually merging files, but need different taxonomy file (since qza is all combined taxon)
#OR troubleshoot qza_to_phylo, follow steps with more control
