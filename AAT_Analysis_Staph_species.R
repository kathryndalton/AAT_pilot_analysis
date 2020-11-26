## Attempt to Get Staph species

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
write_tsv(staph_seq1, path="/Users/kathryndalton/Dropbox/AAT_microbiome/16S/staph_seq_tax", col_names = TRUE)

## Bring Back Manually Added in BLAST Results
staph_seq_tab<-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/staph_seq_tax")

## format into phyloseq uploadable form
taxtable_final<-staph_seq_tab %>% select(-Confidence, -Seq) %>% column_to_rownames("FeatureID") %>% as.matrix()






