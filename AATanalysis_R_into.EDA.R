## AAT Analysis - bring QIIME into R & preliminary EDA

#Prep
setwd("~/Dropbox/AAT_microbiome/16S")
if (!require('tidyverse')) {
  install.packages('tidyverse');
  library(tidyverse);
}
if (!require('devtools')) {
  install.packages('devtools');
  devtools::install_github("jbisanz/qiime2R")
  library(qiime2R);
}
if (!require('qiime2R')) {
  install.packages('devtools');
  devtools::install_github("jbisanz/qiime2R")
  library(qiime2R);
}

#Step 1 upload qza's and metadata
aat_metadata1 <-read_tsv("/Users/kathryndalton/Dropbox/AAT_microbiome/16S/AAT_workbook.txt")
#Import QIIME artifacts - just count table, final ASVs/representative sequences, and taxonomic table
seqtab<-read_qza("aat-table.qza")
names(seqtab)
seqtab$data[1:5,1:5]
class(seqtab)
dim(seqtab)
class(seqtab$data) #seqtab$data is count table not just seqtab
dim(seqtab$data)
taxonomy<-read_qza("aat-taxonomy.qza")
names(taxonomy)
#taxtable<-taxonomy$data %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
#taxtable
class(taxonomy)
dim(taxonomy)
class(taxonomy$data)
dim(taxonomy$data)
taxonomy$data[1:3,1:3]
tree<- read_qza("rooted-tree-aat.qza")
names(tree)
class(tree$data)
dim(tree$data)
#go into specific QIIME diversity metrics @2700rare length - can also upload others eg.5215
setwd("./aat-core-metrics-results-2700/")
shannon<-read_qza("shannon_vector.qza")
head(shannon$data)
UUF_pcoa<-read_qza("unweighted_unifrac_pcoa_results.qza")
UUF_pcoa$data$ProportionExplained[, 1:5]
UUF_pcoa$data$Vectors[1:5, 1:5]
WUF_pcoa<-read_qza("weighted_unifrac_pcoa_results.qza")
FPD<-read_qza("faith_pd_vector.qza")
setwd("~/Dropbox/AAT_microbiome/16S") ##NEED TO GO BACK TO MAIN DIRECTORY

#Merging Metadata w Alpha Data
colnames1 = c("sampleid","shannon")
colnames2 = c("sampleid","fpd")
write.csv(shannon$data, file = "./shannon.csv", sep = " ", row.names = TRUE)
write.csv(FPD$data, file = "./FPD.csv", sep = " ", row.names = TRUE)
shannon_data<-read.csv("./shannon.csv", col.names = colnames1)
FPD_data<-read.csv("./FPD.csv", col.names = colnames2)
m<-merge(aat_metadata1, shannon_data, by="sampleid") #data from only those with shannon/FPD results, not cutoff
data<-merge(m, FPD_data, by="sampleid")
str(data)  #data from all variables even those cutoff
data$interv <- ifelse(data$intervention==0, "Control", 
                             ifelse(data$intervention==1, "Intervention", NA))
#Create Filtered Datasets
human <- filter(data, host=="human")
dog1 <- filter(data, host=="dog")
env <- filter(data, sampletype=="Dust")
control <- filter(data, sampletype=="Extract Control")
blank <- filter(data, sampletype=="Environmental Blank")
nose <- filter(data, sampletype=="Nasopharyngeal swab")

## EDA in R - Alpha numbers
summary(data$shannon)
summary(data$fpd)

by_host <- data %>% group_by(host)
by_sampletype <- data %>% group_by(sampletype)
by_prepost <- data %>% group_by(pre_post)
by_intervention <- data %>% group_by(intervention)
by_samplesite <- data %>% group_by(sample_site)
by_prepost_control <- filter(data, intervention==0) %>% group_by(pre_post)
by_prepost_intervention <- data %>% filter(intervention==1) %>% group_by(pre_post)

by_host_data<-by_host %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                    shannon_med=median(shannon), fpd_med=median(fpd), 
                                    shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                    shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                    shannon_min=min(shannon), fpd_min=min(fpd), 
                                    shannon_max=max(shannon), fpd_max=max(fpd))
by_sampletype_data<-by_sampletype %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                                shannon_med=median(shannon), fpd_med=median(fpd), 
                                                shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                                shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                                shannon_min=min(shannon), fpd_min=min(fpd), 
                                                shannon_max=max(shannon), fpd_max=max(fpd))
by_prepost_data<-by_prepost %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                          shannon_med=median(shannon), fpd_med=median(fpd), 
                                          shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                          shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                          shannon_min=min(shannon), fpd_min=min(fpd), 
                                          shannon_max=max(shannon), fpd_max=max(fpd))
by_intervention_data<-by_intervention %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                                    shannon_med=median(shannon), fpd_med=median(fpd), 
                                                    shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                                    shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                                    shannon_min=min(shannon), fpd_min=min(fpd), 
                                                    shannon_max=max(shannon), fpd_max=max(fpd))
by_samplesite_data<-by_samplesite %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                                shannon_med=median(shannon), fpd_med=median(fpd), 
                                                shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                                shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                                shannon_min=min(shannon), fpd_min=min(fpd), 
                                                shannon_max=max(shannon), fpd_max=max(fpd))
by_prepost_control_data<-by_prepost_control %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                                          shannon_med=median(shannon), fpd_med=median(fpd), 
                                                          shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                                          shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                                          shannon_min=min(shannon), fpd_min=min(fpd), 
                                                          shannon_max=max(shannon), fpd_max=max(fpd))
by_prepost_intervention_data<-by_prepost_intervention %>% summarize(shannon_mean=mean(shannon), fpd_mean=mean(fpd), 
                                                                    shannon_med=median(shannon), fpd_med=median(fpd), 
                                                                    shannon_25=quantile(shannon, 0.25), fpd_25=quantile(fpd, 0.25), 
                                                                    shannon_75=quantile(shannon, 0.75), fpd_75=quantile(fpd, 0.75), 
                                                                    shannon_min=min(shannon), fpd_min=min(fpd), 
                                                                    shannon_max=max(shannon), fpd_max=max(fpd))

list(by_host_data, by_sampletype_data, by_samplesite_data, by_prepost_data, by_intervention_data, by_prepost_control_data, by_prepost_intervention_data)

##Visual EDA alpha
ggplot(data, aes(x=host, y=shannon)) +
  geom_boxplot()
ggplot(data, aes(x=host, y=fpd)) +
  geom_boxplot()
ggplot(data, aes(x=pre_post, y=shannon)) +
  geom_boxplot()
ggplot(human, aes(x=pre_post, y=shannon)) +
  geom_boxplot(aes(fill=interv))
ggplot(human, aes(x=interv, y=shannon)) +
  geom_boxplot(aes(fill=pre_post))
ggplot(dog1, aes(x=interv, y=shannon)) +
  geom_boxplot(aes(fill=pre_post))

ggplot(human, aes(x=pre_post, y=shannon, colour=interv)) +
  geom_point() +
  geom_line(aes(group=subjectid)) + 
  scale_x_discrete(limits=c("PRE", "POST"))

human_prepost<-human %>% group_by(pre_post) %>% summarize(shannon=mean(shannon))
human_prepost_interv<-human %>% group_by(pre_post, interv) %>% summarize(shannon=mean(shannon))
human_prepost$interv <- c("Overall Mean")
colorfills<-c("#F19C4C", "#455E9E", "black")
ggplot(human, aes(x=pre_post, y=shannon, colour=interv)) +
  geom_point(size=1) +
  geom_line(aes(group=subjectid), size=0.5) +
  geom_point(data=human_prepost_interv, size=2) +
  geom_line(data=human_prepost_interv, aes(group=interv), size=1.5) +
  geom_point(data=human_prepost, size=2) +
  geom_line(data=human_prepost, aes(group=interv), size=1.5) +
  scale_colour_manual(values = colorfills) +
  scale_x_discrete(limits=c("PRE", "POST")) +
  ggtitle("Human Alpha Diversity Changes within Visits, by Intervention Group")

ggplot(dog1, aes(x=pre_post, y=shannon, colour=interv)) +
  geom_point(aes(shape=sample_site)) +
  geom_line(aes(group=uniqueid)) + 
  scale_x_discrete(limits=c("PRE", "POST"))

dog_prepost<-dog1 %>% group_by(pre_post) %>% summarize(shannon=mean(shannon))
dog_prepost_interv<-dog1 %>% group_by(pre_post, interv) %>% summarize(shannon=mean(shannon))
dog_prepost$interv <- c("Overall Mean")
colorfills<-c("#F19C4C", "#455E9E", "black")
ggplot(dog1, aes(x=pre_post, y=shannon, colour=interv)) +
  geom_point(aes(shape=sample_site), size=1.5) +
  geom_line(aes(group=uniqueid), size=0.3) +
  geom_point(data=dog_prepost_interv, size=2.5) +
  geom_line(data=dog_prepost_interv, aes(group=interv), size=1.75) +
  geom_point(data=dog_prepost, size=2.5) +
  geom_line(data=dog_prepost, aes(group=interv), size=1.75) +
  scale_colour_manual(values = colorfills) +
  scale_x_discrete(limits=c("PRE", "POST")) +
  ggtitle("Canine Alpha Diversity Changes within Visits, by Intervention Group")

## Alpha differences by contact score and S.a/MRSA status
#Need to bring in new metadata
aat_metadata_child <-read_tsv("AATmetadata1_child_surveyobsculture.txt")
#Change human dataset to wide - one entry per person (pre, post same line)
human_pre <- filter(human, pre_post=="PRE")
human_post <- filter(human, pre_post=="POST")
human1<-subset(human_pre, select = c(subjectid, shannon, fpd, batch))
human2<-subset(human_post, select = c(subjectid, shannon, fpd, batch))
human3<-subset(human, select = c(sampletype, subjectid, host, record_ID, 
                                 visit_number, match_week, sample_site, intervention, interv))
human4<-unique(human3)
human5<-merge(human4, human1, by="subjectid")
human6<-merge(human5, human2, by="subjectid", suffixes = c("_pre", "_post"))
#Merge just subsetted human dataset to other metadata
child<-merge(human6, aat_metadata_child, by="subjectid")
#generate alpha change variable
child$shannon_change = child$shannon_post - child$shannon_pre
hist(child$shannon_change)
summary(child$shannon_change)
group_by(child, interv) %>% summarize(shannon_change=mean(shannon_change))
group_by(child, contactscore) %>% summarize(shannon_change=mean(shannon_change))
ggplot(child, aes(x=interv, y=shannon_change)) +
  geom_boxplot()
ggplot(child, aes(x=contactscore, y=shannon_change)) +
  geom_boxplot()
ggplot(child, aes(x=contactscore, y=shannon_change)) +
  geom_boxplot(aes(fill=interv))
child$fpd_change = child$fpd_post - child$fpd_pre
hist(child$fpd_change)
summary(child$fpd_change)
group_by(child, interv) %>% summarize(fpd_change=mean(fpd_change))
group_by(child, contactscore) %>% summarize(fpd_change=mean(fpd_change))
ggplot(child, aes(x=interv, y=fpd_change)) +
  geom_boxplot()
ggplot(child, aes(x=contactscore, y=fpd_change)) +
  geom_boxplot()
ggplot(child, aes(x=contactscore, y=fpd_change)) +
  geom_boxplot(aes(fill=interv))
wilcox.test(shannon_change ~ contactscore, data = child, paired = FALSE) #Mann-Whitney U Test
wilcox.test(shannon_change ~ interv, data = child, paired = FALSE)
wilcox.test(fpd_change ~ contactscore, data = child, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = child, paired = FALSE)
child_control <- filter(child, interv=="Control")
child_interv<-filter(child, interv=="Intervention")
child_highcontact<-filter(child, contactscore=="High Contact")
child_lowcontact<-filter(child, contactscore=="Low contact")
wilcox.test(shannon_change ~ contactscore, data = child_control, paired = FALSE)
wilcox.test(shannon_change ~ contactscore, data = child_interv, paired = FALSE)
wilcox.test(fpd_change ~ contactscore, data = child_control, paired = FALSE)
wilcox.test(fpd_change ~ contactscore, data = child_interv, paired = FALSE)
wilcox.test(shannon_change ~ interv, data = child_highcontact, paired = FALSE)
wilcox.test(shannon_change ~ interv, data = child_lowcontact, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = child_highcontact, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = child_lowcontact, paired = FALSE)

## Shannon Alpha Difference by Culture Results
group_by(child, kidsa_pre) %>% summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
group_by(child, kidmrsa_pre) %>% summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
group_by(child, kidsa_post) %>% summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))
group_by(child, kidmrsa_post) %>% summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))
ggplot(child, aes(x=kidsa_pre, y=shannon_pre)) +
  geom_boxplot(aes(group=kidsa_pre))
ggplot(child, aes(x=kidmrsa_pre, y=shannon_pre)) +
  geom_boxplot(aes(group=kidmrsa_pre))
ggplot(child, aes(x=kidsa_post, y=shannon_post)) +
  geom_boxplot(aes(group=kidsa_post))
ggplot(child, aes(x=kidmrsa_post, y=shannon_post)) +
  geom_boxplot(aes(group=kidmrsa_post))
child$kidsa_gain<-ifelse(child$kidsa_pre==0&child$kidsa_post==1, 1, 0)
child$kidmrsa_gain<-ifelse(child$kidmrsa_pre==0&child$kidmrsa_post==1, 1, 0)
ggplot(child, aes(x=kidsa_gain, y=shannon_change)) +
  geom_boxplot(aes(group=kidsa_gain)) +
  geom_point(aes(color=interv))
ggplot(child, aes(x=kidmrsa_gain, y=shannon_change)) +
  geom_boxplot(aes(group=kidmrsa_gain)) +
  geom_point(aes(color=interv))
ggplot(child, aes(x=kidsa_gain, y=shannon_change)) +
  geom_boxplot(aes(group=kidsa_gain)) +
  geom_point(aes(color=contactscore))
ggplot(child, aes(x=kidmrsa_gain, y=shannon_change)) +
  geom_boxplot(aes(group=kidmrsa_gain)) +
  geom_point(aes(color=contactscore))
wilcox.test(shannon_change ~ kidsa_gain, data = child, paired = FALSE)
wilcox.test(shannon_change ~ kidmrsa_gain, data = child, paired = FALSE)
child_control <- filter(child, interv=="Control")
child_interv<-filter(child, interv=="Intervention")
child_highcontact<-filter(child, contactscore=="High Contact")
child_lowcontact<-filter(child, contactscore=="Low contact")
wilcox.test(shannon_change ~ kidsa_gain, data = child_control, paired = FALSE)
wilcox.test(shannon_change ~ kidsa_gain, data = child_interv, paired = FALSE)
wilcox.test(shannon_change ~ kidmrsa_gain, data = child_control, paired = FALSE)
wilcox.test(shannon_change ~ kidmrsa_gain, data = child_interv, paired = FALSE) ##doesnt work because no mrs neg in interveration
wilcox.test(shannon_change ~ kidsa_gain, data = child_highcontact, paired = FALSE)
wilcox.test(shannon_change ~ kidsa_gain, data = child_lowcontact, paired = FALSE)
wilcox.test(shannon_change ~ kidmrsa_gain, data = child_highcontact, paired = FALSE)
wilcox.test(shannon_change ~ kidmrsa_gain, data = child_lowcontact, paired = FALSE)

## Fpd Alpha Difference by Culture Results
group_by(child, kidsa_pre) %>% summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
group_by(child, kidmrsa_pre) %>% summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
group_by(child, kidsa_post) %>% summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))
group_by(child, kidmrsa_post) %>% summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))
ggplot(child, aes(x=kidsa_pre, y=fpd_pre)) +
  geom_boxplot(aes(group=kidsa_pre))
ggplot(child, aes(x=kidmrsa_pre, y=fpd_pre)) +
  geom_boxplot(aes(group=kidmrsa_pre))
ggplot(child, aes(x=kidsa_post, y=fpd_post)) +
  geom_boxplot(aes(group=kidsa_post))
ggplot(child, aes(x=kidmrsa_post, y=fpd_post)) +
  geom_boxplot(aes(group=kidmrsa_post))
child$kidsa_gain<-ifelse(child$kidsa_pre==0&child$kidsa_post==1, 1, 0)
child$kidmrsa_gain<-ifelse(child$kidmrsa_pre==0&child$kidmrsa_post==1, 1, 0)
ggplot(child, aes(x=kidsa_gain, y=fpd_change)) +
  geom_boxplot(aes(group=kidsa_gain)) +
  geom_point(aes(color=interv))
ggplot(child, aes(x=kidmrsa_gain, y=fpd_change)) +
  geom_boxplot(aes(group=kidmrsa_gain)) +
  geom_point(aes(color=interv))
ggplot(child, aes(x=kidsa_gain, y=fpd_change)) +
  geom_boxplot(aes(group=kidsa_gain)) +
  geom_point(aes(color=contactscore))
ggplot(child, aes(x=kidmrsa_gain, y=fpd_change)) +
  geom_boxplot(aes(group=kidmrsa_gain)) +
  geom_point(aes(color=contactscore))
wilcox.test(fpd_change ~ kidsa_gain, data = child, paired = FALSE)
wilcox.test(fpd_change ~ kidmrsa_gain, data = child, paired = FALSE)
child_control <- filter(child, interv=="Control")
child_interv<-filter(child, interv=="Intervention")
child_highcontact<-filter(child, contactscore=="High Contact")
child_lowcontact<-filter(child, contactscore=="Low contact")
wilcox.test(fpd_change ~ kidsa_gain, data = child_control, paired = FALSE)
wilcox.test(fpd_change ~ kidsa_gain, data = child_interv, paired = FALSE)
wilcox.test(fpd_change ~ kidmrsa_gain, data = child_control, paired = FALSE)
wilcox.test(fpd_change ~ kidmrsa_gain, data = child_interv, paired = FALSE) ##doesnt work because no mrs neg in interveration
wilcox.test(fpd_change ~ kidsa_gain, data = child_highcontact, paired = FALSE)
wilcox.test(fpd_change ~ kidsa_gain, data = child_lowcontact, paired = FALSE)
wilcox.test(fpd_change ~ kidmrsa_gain, data = child_highcontact, paired = FALSE)
wilcox.test(fpd_change ~ kidmrsa_gain, data = child_lowcontact, paired = FALSE)


## Alpha differences by culture in Dogs

#Change dog dataset to wide - one entry per dog SITE N, MO ... (pre, post same line)
dog_pre <- filter(dog, pre_post=="PRE")
dog_post <- filter(dog, pre_post=="POST")
dog1<-subset(dog_pre, select = c(uniqueid, shannon, fpd, batch))
dog2<-subset(dog_post, select = c(uniqueid, shannon, fpd, batch))
dog3<-subset(dog, select = c(sampletype, subjectid, host, record_ID, 
                             visit_number, match_week, sample_site, uniqueid, intervention, interv))
dog4<-unique(dog3)
dog5<-merge(dog4, dog1, by="uniqueid")
dog6<-merge(dog5, dog2, by="uniqueid", suffixes = c("_pre", "_post"))
#Merge just subsetted human dataset to other metadata
dog_culture <-read_tsv("AATmetadata2_dog_culture.txt")
dogw<-merge(dog6, dog_culture, by="match_week")

#generate alpha change variable
dogw$shannon_change <- dog$shannon_post - dog$shannon_pre
hist(dogw$shannon_change)
summary(dogw$shannon_change)
group_by(dogw, interv) %>% summarize(shannon_change=mean(shannon_change))
group_by(dogw, sampletype) %>% summarize(shannon_change=mean(shannon_change))
ggplot(dogw, aes(x=interv, y=shannon_change)) +
  geom_boxplot()
ggplot(dogw, aes(x=sample_site, y=shannon_change)) +
  geom_boxplot()
ggplot(dogw, aes(x=sample_site, y=shannon_change)) +
  geom_boxplot(aes(fill=interv))

dogw$fpd_change = dog$fpd_post - dog$fpd_pre
hist(dogw$fpd_change)
summary(dogw$fpd_change)
group_by(dogw, interv) %>% summarize(fpd_change=mean(fpd_change))
group_by(dogw, sample_site) %>% summarize(fpd_change=mean(fpd_change))
ggplot(dogw, aes(x=interv, y=fpd_change)) +
  geom_boxplot()
ggplot(dogw, aes(x=sample_site, y=fpd_change)) +
  geom_boxplot()
ggplot(dogw, aes(x=sample_site, y=fpd_change)) +
  geom_boxplot(aes(fill=interv))

kruskal.test(shannon_change ~ sample_site, data = dogw) #Kruskal-Wallis RankSum Test >2categories
wilcox.test(shannon_change ~ interv, data = dogw, paired = FALSE) #Mann-Whitney U Test
kruskal.test(fpd_change ~ sample_site, data = dogw)
wilcox.test(fpd_change ~ interv, data = dogw, paired = FALSE)
dog_control <- filter(dogw, interv=="Control")
dog_interv<-filter(dogw, interv=="Intervention")
dog_nasal<-filter(dogw, sample_site=="N")
dog_oral<-filter(dogw, sample_site=="MO")
dog_inguinal<-filter(dogw, sample_site=="I")
dog_perineal<-filter(dogw, sample_site=="R")
kruskal.test(shannon_change ~ sample_site, data = dog_control)
kruskal.test(shannon_change ~ sample_site, data = dog_interv)
kruskal.test(fpd_change ~ sample_site, data = dog_control)
kruskal.test(fpd_change ~ sample_site, data = dog_interv)
wilcox.test(shannon_change ~ interv, data = dog_nasal, paired = FALSE)
wilcox.test(shannon_change ~ interv, data = dog_oral, paired = FALSE)
wilcox.test(shannon_change ~ interv, data = dog_inguinal, paired = FALSE)
wilcox.test(shannon_change ~ interv, data = dog_perineal, paired = FALSE) ##p<0.05!
wilcox.test(fpd_change ~ interv, data = dog_nasal, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = dog_oral, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = dog_inguinal, paired = FALSE)
wilcox.test(fpd_change ~ interv, data = dog_perineal, paired = FALSE) ##p<0.05!

#Site specific - Shannon
group_by(dogw, doganysa_pre, doganymrsa_pre) %>% summarize(shannon_pre_mean=mean(shannon_pre), 
                                                           shannon_pre_median=median(shannon_pre))
filter(dogw, sample_site=="N") %>% group_by(dognsa_pre, dognmrsa_pre) %>% 
  summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
filter(dogw, sample_site=="MO") %>% group_by(dogmosa_pre, dogmomrsa_pre) %>% 
  summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
filter(dogw, sample_site=="I") %>% group_by(dogisa_pre, dogimrsa_pre) %>% 
  summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
filter(dogw, sample_site=="R") %>% group_by(dogrsa_pre, dogrmrsa_pre) %>% 
  summarize(shannon_pre_mean=mean(shannon_pre), shannon_pre_median=median(shannon_pre))
group_by(dogw, doganysa_post, doganymrsa_post) %>% summarize(shannon_post_mean=mean(shannon_post), 
                                                             shannon_post_median=median(shannon_post))
filter(dogw, sample_site=="N") %>% group_by(dognsa_post, dognmrsa_post) %>% 
  summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))
filter(dogw, sample_site=="MO") %>% group_by(dogmosa_post, dogmomrsa_post) %>% 
  summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))
filter(dogw, sample_site=="I") %>% group_by(dogisa_post, dogimrsa_post) %>% 
  summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))
filter(dogw, sample_site=="R") %>% group_by(dogrsa_post, dogrmrsa_post) %>% 
  summarize(shannon_post_mean=mean(shannon_post), shannon_post_median=median(shannon_post))

ggplot(dogw, aes(x=shannon_pre, y=shannon_post)) + 
  geom_point(aes(colour=sample_site, shape=sample_site))
ggplot(dogw, aes(x=shannon_pre, y=shannon_post)) + 
  geom_point(aes(colour=doganysa_pre, shape=sample_site))
ggplot(dogw, aes(x=shannon_pre, y=shannon_post)) + 
  geom_point(aes(colour=doganymrsa_pre, shape=sample_site))

dogw$dogany_gain<-ifelse(dog$doganysa_pre==0&dog$doganysa_post==1, "S.a. Gain", 
                         ifelse(dog$doganymrsa_pre==0&dog$doganymrsa_post==1, "MRSA Gain", "No"))
dogw$dogn_gain<-ifelse(dog$dognsa_pre==0&dog$dognsa_post==1, "S.a. Gain", 
                       ifelse(dog$dognmrsa_pre==0&dog$dognmrsa_post==1, "MRSA Gain", "No"))
dogw$dogmo_gain<-ifelse(dog$dogmosa_pre==0&dog$dogmosa_post==1, "S.a. Gain", 
                        ifelse(dog$dogmomrsa_pre==0&dog$dogmomrsa_post==1, "MRSA Gain", "No"))
dogw$dogi_gain<-ifelse(dog$dogisa_pre==0&dog$dogisa_post==1, "S.a. Gain", 
                       ifelse(dog$dogimrsa_pre==0&dog$dogimrsa_post==1, "MRSA Gain", "No"))
dogw$dogr_gain<-ifelse(dog$dogrsa_pre==0&dog$dogrsa_post==1, "S.a. Gain", 
                       ifelse(dog$dogrmrsa_pre==0&dog$dogrmrsa_post==1, "MRSA Gain", "No"))
a.s<-ggplot(dogw, aes(x=dogany_gain, y=shannon_change)) +
  geom_boxplot() +
  ggtitle("All")  
n.s<-filter(dogw, sample_site=="N") %>% ggplot(aes(x=dogn_gain, y=shannon_change)) + 
  geom_boxplot() +
  ggtitle("Nose")
mo.s<-filter(dogw, sample_site=="MO") %>% ggplot(aes(x=dogmo_gain, y=shannon_change)) + 
  geom_boxplot() +
  ggtitle("Oral")
i.s<-filter(dogw, sample_site=="I") %>% ggplot(aes(x=dogi_gain, y=shannon_change)) + 
  geom_boxplot() +
  ggtitle("Inguinal")    
r.s<-filter(dogw, sample_site=="R") %>% ggplot(aes(x=dogr_gain, y=shannon_change)) + 
  geom_boxplot() +
  ggtitle("Perineal")
gridExtra::grid.arrange(a.s, n.s, mo.s, i.s, r.s, ncol=2, nrow=3)

#Site specific - Fpd
group_by(dogw, doganysa_pre, doganymrsa_pre) %>% summarize(fpd_pre_mean=mean(fpd_pre), 
                                                           fpd_pre_median=median(fpd_pre))
filter(dogw, sample_site=="N") %>% group_by(dognsa_pre, dognmrsa_pre) %>% 
  summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
filter(dogw, sample_site=="MO") %>% group_by(dogmosa_pre, dogmomrsa_pre) %>% 
  summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
filter(dogw, sample_site=="I") %>% group_by(dogisa_pre, dogimrsa_pre) %>% 
  summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
filter(dogw, sample_site=="R") %>% group_by(dogrsa_pre, dogrmrsa_pre) %>% 
  summarize(fpd_pre_mean=mean(fpd_pre), fpd_pre_median=median(fpd_pre))
group_by(dogw, doganysa_post, doganymrsa_post) %>% summarize(fpd_post_mean=mean(fpd_post), 
                                                             fpd_post_median=median(fpd_post))
filter(dogw, sample_site=="N") %>% group_by(dognsa_post, dognmrsa_post) %>% 
  summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))
filter(dogw, sample_site=="MO") %>% group_by(dogmosa_post, dogmomrsa_post) %>% 
  summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))
filter(dogw, sample_site=="I") %>% group_by(dogisa_post, dogimrsa_post) %>% 
  summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))
filter(dogw, sample_site=="R") %>% group_by(dogrsa_post, dogrmrsa_post) %>% 
  summarize(fpd_post_mean=mean(fpd_post), fpd_post_median=median(fpd_post))

ggplot(dogw, aes(x=fpd_pre, y=fpd_post)) + 
  geom_point(aes(colour=sample_site, shape=sample_site))
ggplot(dogw, aes(x=fpd_pre, y=fpd_post)) + 
  geom_point(aes(colour=doganysa_pre, shape=sample_site))
ggplot(dogw, aes(x=fpd_pre, y=fpd_post)) + 
  geom_point(aes(colour=doganymrsa_pre, shape=sample_site))

a.f<-ggplot(dogw, aes(x=dogany_gain, y=fpd_change)) +
  geom_boxplot() +
  ggtitle("All")  
n.f<-filter(dogw, sample_site=="N") %>% ggplot(aes(x=dogn_gain, y=fpd_change)) + 
  geom_boxplot() +
  ggtitle("Nose")
mo.f<-filter(dogw, sample_site=="MO") %>% ggplot(aes(x=dogmo_gain, y=fpd_change)) + 
  geom_boxplot() +
  ggtitle("Oral")
i.f<-filter(dogw, sample_site=="I") %>% ggplot(aes(x=dogi_gain, y=fpd_change)) + 
  geom_boxplot() +
  ggtitle("Inguinal")    
r.f<-filter(dogw, sample_site=="R") %>% ggplot(aes(x=dogr_gain, y=fpd_change)) + 
  geom_boxplot() +
  ggtitle("Perineal")
gridExtra::grid.arrange(a.f, n.f, mo.f, i.f, r.f, ncol=2, nrow=3)


## Beta Diversity
UUF_pcoa$data$ProportionExplained[, 1:5]
UUF_pcoa$data$Vectors[1:5, 1:5]

#Step 1 - format table
Data_UUF_pcoa<-merge(data, UUF_pcoa$data$Vectors, by.x="sampleid", by.y = "SampleID")
Data_WUF_pcoa<-merge(data, WUF_pcoa$data$Vectors, by.x="sampleid", by.y = "SampleID")
Data_UUF_pcoa_kid<-filter(Data_UUF_pcoa, host=="human")
Data_WUF_pcoa_kid<-filter(Data_WUF_pcoa, host=="human")
Data_UUF_pcoa_dog<-filter(Data_UUF_pcoa, host=="dog")
Data_WUF_pcoa_dog<-filter(Data_WUF_pcoa, host=="dog")
#Plots
ggplot(Data_UUF_pcoa_kid, aes(x=PC1, y=PC2, shape=pre_post)) +
  geom_point(aes(color=interv), size=2) +
  geom_line(aes(group=subjectid, color=interv)) +
  ggtitle("Kid Unweighted Unifrac PCOA")
ggplot(Data_WUF_pcoa_kid, aes(x=PC1, y=PC2, shape=pre_post)) +
  geom_point(aes(color=interv), size=2) +
  geom_line(aes(group=subjectid, color=interv)) +
  ggtitle("Kid Weighted Unifrac PCOA")
ggplot(Data_UUF_pcoa_dog, aes(x=PC1, y=PC2, shape=sample_site)) +
  geom_point(aes(color=interv), size=2) +
  geom_line(aes(group=uniqueid, color=interv)) +
  ggtitle("Dog Unweighted Unifrac PCOA")
ggplot(Data_WUF_pcoa_dog, aes(x=PC1, y=PC2, shape=sample_site)) +
  geom_point(aes(color=interv), size=2) +
  geom_line(aes(group=uniqueid, color=interv)) +
  ggtitle("Dog Weighted Unifrac PCOA")
ggplot(Data_UUF_pcoa_kid, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=shannon, color=interv)) +
  ggtitle("Kid Unweighted Unifrac PCOA with Alpha")
ggplot(Data_WUF_pcoa_kid, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=shannon, color=interv)) +
  ggtitle("Dog Weighted Unifrac PCOA with Alpha")
ggplot(Data_UUF_pcoa_dog, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=shannon, color=interv, shape=sample_site)) +
  ggtitle("Kid Unweighted Unifrac PCOA with Alpha")
ggplot(Data_WUF_pcoa_dog, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=shannon, color=interv, shape=sample_site)) +
  ggtitle("Dog Weighted Unifrac PCOA with Alpha")
filter(Data_UUF_pcoa, !is.na(match_week)) %>% ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=factor(match_week), shape=host)) +
  scale_color_discrete() + labs(color="Visit Week") +
  ggtitle("Uniweighted Unifrace by Visit")
filter(Data_WUF_pcoa, !is.na(match_week)) %>% ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=factor(match_week), shape=host)) +
  scale_color_discrete() + labs(color="Visit Week") +
  ggtitle("Weighted Unifrace by Visit")

##Beta Differences
uuf_distance <-read_tsv("uuf-distance.tsv")
wuf_distance <-read_tsv("wuf-distance.tsv")
u1<-merge(data, uuf_distance, by.x="sampleid", by.y="SubjectID1", all=TRUE)
u2<-merge(u1, data, by.x="SubjectID2", by.y="sampleid", all=TRUE) 
#sampleid becomes Subject1 =x and Subject2 =y
w1<-merge(data, wuf_distance, by.x="sampleid", by.y="SubjectID1", all=TRUE)
w2<-merge(w1, data, by.x="SubjectID2", by.y="sampleid", all=TRUE) 
#sampleid becomes Subject1 =x and Subject2 =y
u3<-filter(u2, !(is.na(interv.x)==TRUE)) 
u4<-filter(u3, !(is.na(interv.y)==TRUE)) #removes controls
w3<-filter(w2, !(is.na(interv.x)==TRUE))
w4<-filter(w3, !(is.na(interv.y)==TRUE)) #NAs are ENV dust
u4$host_match<-as.factor(paste(u4$host.x, u4$host.y, sep = "-"))
w4$host_match<-as.factor(paste(w4$host.x, w4$host.y, sep = "-"))
u4$hostmatch<-ifelse(u4$host_match=="human-human", "kid-kid", 
                     ifelse(u4$host_match=="human-dog", "kid-dog",
                            ifelse(u4$host_match=="dog-human", "kid-dog",
                                   ifelse(u4$host_match=="human-NA", "kid-env",
                                          ifelse(u4$host_match=="NA-human", "kid-env",
                                                 ifelse(u4$host_match=="dog-NA", "dog-env",
                                                        ifelse(u4$host_match=="NA-dog", "dog-env", 
                                                               ifelse(u4$host_match=="dog-dog", "dog-dog",
                                                                      ifelse(u4$host_match=="NA-NA", "env-env", NA)))))))))  
w4$hostmatch<-ifelse(w4$host_match=="human-human", "kid-kid", 
                     ifelse(w4$host_match=="human-dog", "kid-dog",
                            ifelse(w4$host_match=="dog-human", "kid-dog",
                                   ifelse(w4$host_match=="human-NA", "kid-env",
                                          ifelse(w4$host_match=="NA-human", "kid-env",
                                                 ifelse(w4$host_match=="dog-NA", "dog-env",
                                                        ifelse(w4$host_match=="NA-dog", "dog-env", 
                                                               ifelse(w4$host_match=="dog-dog", "dog-dog",
                                                                      ifelse(w4$host_match=="NA-NA", "env-env", NA)))))))))  
uuf_matchweek<-subset(u4, match_week.x==match_week.y)
wuf_matchweek<-subset(w4, match_week.x==match_week.y)
filter(uuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(hostmatch) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
filter(wuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(hostmatch) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
filter(uuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  kruskal.test(Distance ~ hostmatch)
filter(uuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  kruskal.test(Distance ~ visit_number.x)
filter(wuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  kruskal.test(Distance ~ hostmatch)
filter(wuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  kruskal.test(Distance ~ visit_number.x)

filter(uuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  ggplot(aes(x=visit_number.x, y=Distance, fill=hostmatch)) + geom_boxplot() +ggtitle("Unweighted Unifrac")
filter(wuf_matchweek, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  ggplot(aes(x=visit_number.x, y=Distance, fill=hostmatch)) + geom_boxplot() +ggtitle("Weighted Unifrac")

u4$visit_match<-ifelse(u4$match_week.x==u4$match_week.y, "Same Visit", "Different Visit")
w4$visit_match<-ifelse(w4$match_week.x==w4$match_week.y, "Same Visit", "Different Visit")
filter(u4, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  ggplot(aes(x=hostmatch, y=Distance, fill=visit_match)) + geom_boxplot() +ggtitle("Unweighted Unifrac")
filter(w4, pre_post.x=="PRE" & pre_post.y=="POST") %>% 
  ggplot(aes(x=hostmatch, y=Distance, fill=visit_match)) + geom_boxplot() +ggtitle("Weighted Unifrac")

#Beta Differences Human-Dog by Contact Score

uuf_human<-filter(u4, host.x=="human")
wuf_human<-filter(w4, host.x=="human")
uh<-merge(uuf_human, child, by.x ="subjectid.x", by.y="subjectid")
wh<-merge(wuf_human, child, by.x ="subjectid.x", by.y="subjectid")
uh1<-merge(Data_UUF_pcoa_kid, child, by="subjectid")
wh1<-merge(Data_WUF_pcoa_kid, child, by="subjectid")
filter(uh, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(visit_match) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
filter(uh, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(match_week.x.x) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
filter(wh, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(visit_match) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
filter(wh, pre_post.x=="PRE" & pre_post.y=="POST") %>%
  group_by(match_week.x.x) %>% 
  summarize(dist_mean=mean(Distance), dist_median=median(Distance))
ggplot(uh1, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour=factor(match_week), shape=contactscore)) +
  scale_color_discrete() +
  ggtitle("Kid Unweighted Unifrac PCOA, Contact Score and Visit")
ggplot(uh1, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour=factor(dog), shape=contactscore)) +
  scale_color_discrete() + 
  ggtitle("Kid Unweighted Unifrac PCOA, Contact Score and Dog")
ggplot(wh1, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour=factor(match_week), shape=contactscore)) +
  scale_color_discrete() +
  ggtitle("Kid Weighted Unifrac PCOA, Contact Score and Visit")
ggplot(wh1, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour=factor(dog), shape=contactscore)) +
  scale_color_discrete() + 
  ggtitle("Kid Weighted Unifrac PCOA, Contact Score and Dog")


