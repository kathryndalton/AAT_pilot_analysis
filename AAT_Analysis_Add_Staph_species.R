## Attempt to Get Staph species
## Adding in S. pseud and S. schliefieri to greengenes databased by using BLAST ##

##Set Up
seqs<-staph_tab<-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/sequences.tsv")
staph_tab<-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/staph_tab.tsv")
staph_taxtable1<-staph_tab %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
staph_taxtable2<-subset(staph_taxtable1, staph_taxtable1$Genus=="g__Staphylococcus")
staph_taxtable<-subset(staph_taxtable2, staph_taxtable2$Species =="s__")
seq1<-rownames_to_column(seqs)
is.odd <- function(v) v %% 2 != 0  
seq1<-transform(seq1, rowname = as.numeric(rowname))
class(seq1$rowname)
seq1$odd<-is.odd(seq1$rowname)
seq1.name<-subset(seq1, odd==TRUE)
seq1.seq<-subset(seq1, odd==FALSE)
seq1.name<-transmute(seq1.name, FeatureID=Seq)
seq1.name<-rownames_to_column(seq1.name)
seq1.seq$row<-c(1:14183)
seq2<-merge(x=seq1.seq, y=seq1.name, by.x="row", by.y="rowname", all=TRUE)

#Isolate Staph Sequences
staph_seq<-merge(x=seq2, y=staph_taxtable, by.x="FeatureID", by.y="Feature ID", all.y=TRUE, all.x=FALSE)
staph_seq1<-subset(staph_seq, is.na(staph_seq$Seq)==FALSE)
staph_seq1$row<-NULL
staph_seq1$rowname<-NULL
staph_seq1$odd<-NULL
staph_seq2<-subset(staph_seq1, select=Seq:FeatureID)
staph_seq2$start<-c(">")
staph_seq2<-staph_seq2[, c("start", "FeatureID", "Seq")]
write_tsv(staph_seq1, path="/Users/kathryndalton/Dropbox/AAT_microbiome/16S/staph_seq_tax_pre", col_names = TRUE)

#BLAST results individually - for future need to figure out to get into format to upload as one file

## Bring Back Manually Added in BLAST Results
staph_seq_tab<-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/staph_seq_tax")
staph_seq_tab$Seq<-NULL
#need to merge w/ staph_taxtable1 - removing suplicat entries
staph1<-full_join(x=staph_seq_tab, y=staph_taxtable1, by =c("FeatureID"="Feature ID"), match="first")

staph1$Kingdom<-ifelse(is.na(staph1$Kingdom.x), staph1$Kingdom.y, staph1$Kingdom.x)
staph1$Phylum<-ifelse(is.na(staph1$Phylum.x), staph1$Phylum.y, staph1$Phylum.x)
staph1$Class<-ifelse(is.na(staph1$Class.x), staph1$Class.y, staph1$Class.x)
staph1$Order<-ifelse(is.na(staph1$Order.x), staph1$Order.y, staph1$Order.x)
staph1$Family<-ifelse(is.na(staph1$Family.x), staph1$Family.y, staph1$Family.x)
staph1$Genus<-ifelse(is.na(staph1$Genus.x), staph1$Genus.y, staph1$Genus.x)
staph1$Species<-ifelse(is.na(staph1$Species.x), staph1$Species.y, staph1$Species.x)
staph1$Confidence<-ifelse(is.na(staph1$Confidence.x), staph1$Confidence.y, staph1$Confidence.x)

taxtable_final<-subset(staph1, select= -c(Kingdom.x:Confidence.y))
write_tsv(as.data.frame(taxtable_final), "/Users/kathryndalton/Dropbox/AAT_microbiome/16S/aat-taxtable_final.tsv") 

#Bring back to phyloseq object
seqtab<-read_qza("aat-featuretable-final.qza")
tree<-read_qza("rooted-tree-aat.qza")
aat_metadata <-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/AAT_workbook_contact.txt") 
phylo.final<-phyloseq(
  otu_table(seqtab$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable_final) %>% select(-Confidence) %>% column_to_rownames("FeatureID") %>% as.matrix()),
  sample_data(aat_metadata %>% as.data.frame() %>% column_to_rownames("SampleID"))
)
phylo.final
