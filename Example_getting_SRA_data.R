# Getting PETS data from NCBI SRA website into fastq form#

#SRA Study = SRP042152, 240 runs
#single PCR amplicon, Illumina MiSeq#
SRR1548260
SRX687085

#Install BiocManager and SRAdb# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  +     install.packages("BiocManager")
BiocManager::install("SRAdb")
library(SRAdb)

##dowload SQLite database to interact with SRAdb
sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
str(sra_dbname)
str(sra_con)
# OR
sqlfile = getSRAdbFile()
sra_con1 = dbConnect(SQLite(), sqlfile)

## Get column descriptions
a <- colDescriptions(sra_con=sra_con)[1:5,]
head(a)
View(a)
sra_tables <- dbListTables(sra_con)
sra_tables
dbListFields(sra_con,"study")
colDesc <- colDescriptions(sra_con=sra_con)[1:5,]
colDesc[, 1:4]

## List fastq file ftp or fasp addresses associated with "SRP042152"
listSRAfile (in_acc = c("SRP042152"), sra_con = sra_con, fileType = 'sra') 
listSRAfile (in_acc = c("SRP042152"), sra_con = sra_con, fileType = 'sra', srcType='fasp')

PETS_fastq_list1 = listSRAfile( c("SRP042152"), sra_con, fileType = 'sra' )
PETS_fastq_info1 = getSRAinfo( c("SRP042152"), sra_con, sraType = 'sra' )

listSRAfile ( c("SRP042152"), sra_con, fileType ='fastq', srcType='fasp')
getFASTQinfo( c("SRP042152"), sra_con, srcType = 'fasp' )


#list fastq file names including ftp addresses
#listSRAfile function does not check file availability, size and date of the fastq files on the server
PETS_fastq_list <- listSRAfile("SRP042152", sra_con = sra_con, fileType = "fastq",
                  srcType = "ftp")
nrow(PETS_fastq_list)

#commands get fastq file information and essential associated meta data from EBI ENA web site
PETS_fastq_info1 <- getFASTQinfo(in_acc = c("SRX000122"), srcType = "ftp")

## Get file size and date from NCBI site for available fastq files associated with "SRP042152" #
getSRAinfo (in_acc=c("SRP042152"), sra_con=sra_con, sraType='sra')

PETS_fastq_info2 <- dbGetQuery(sra_con, paste("SELECT library strategy AS ’Library Strategy’,",
                                "count( * ) AS Runs FROM ‘experiment‘", "GROUP BY library strategy order by Runs DESC",
                                sep = " "))
head(rs)

##Getting Data
PETS_fastq_files <- getSRAfile('SRP042152',sra_con,fileType='sra')
SRR1548663
PETS_fastq_files <- getSRAfile('SRR1548663',sra_con,fileType='sra')

getSRAfile( in_acc = c("SRP042152"), 
            sra_con = sra_con,
            destDir = getwd(), 
            fileType = 'sra', 
            srcType = 'ftp', 
            makeDirectory = FALSE,
            method = 'curl',
            ascpCMD = NULL)
