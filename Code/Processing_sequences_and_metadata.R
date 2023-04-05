#Code to process DENV 2 from sequence and metadata downloaded from GenBank for Vietnam and SEA 
###Summary of Code###
#Processes metadata 
#Selects for relevant countries from SEA
#Matches metadata to fasta sequences for naming

###Code Start###

setwd("/Users/rhysinward/Documents/Vietnam_nextstrain/Cosmo_manuscript/updated_run")

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main functions to run script
files.sources = list.files(path = "./Main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./Main/", files.sources[i]), collapse = ''))
}

#Load packages

library(tidyverse)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr) #may need to load development version
library(stringr)   
library(tidyr) 
library(seqinr)
library(ggplot2)
library(ape)
library(purrr)
library(stringr)

###DENV-2###

#Need to fill in dates with yyy-mm-dd
#need to split GEO_location to obtain states
#need to use length to obtain e gene vs full length 
#need to select for only SEA countries 

#load in metadata

metadata_DENV2 = as.data.frame(fread('Data/sequences.csv',drop = c('Submitters','Species',
                                                                                  'Molecule_type','Sequence_Type',
                                                                                  'Genotype','Segment','USA',
                                                                                  'Isolation_Source','Release_Date',
                                                                   'Organization', 'Org_location','Org_location')))

#Process dates + check dates

metadata_DENV2 <- process_date(metadata_DENV2)

#select for relevant countries + host

metadata_DENV2$Country[metadata_DENV2$Country == "Timor-Leste"] <- "East Timor"

Relevant_countries <- c("Viet Nam", "Sri Lanka", "China", "Indonesia", "Malaysia",
                        "Myanmar", "Cambodia", "India", "Philippines", "Thailand",
                        "Laos","Borneo", "Bangladesh", "Bhutan", "Singapore", "Taiwan",
                        "East Timor","Nepal")

metadata_DENV2 <- metadata_DENV2 %>%
  filter(Country %in% Relevant_countries)  %>%
  filter(Host == 'Homo sapiens') %>%
  filter(Date >= '2010-01-01')

#extract state level information

metadata_DENV2 <- metadata_DENV2 %>% separate(Geo_Location, c("Country","State"), sep = ":")
metadata_DENV2$State <- trimws(metadata_DENV2$State)

#Classifying dengue genomes as WG/PG/E gene

metadata_DENV2$Sequence.Type <- NA

for (i in 1:nrow(metadata_DENV2)) {
  if(metadata_DENV2$Length[i] <= 3000){
    metadata_DENV2$Sequence.Type[i] <- "E"
  } else if (metadata_DENV2$Length[i] > 3000 & metadata_DENV2$Length[i] <= 8500){
    metadata_DENV2$Sequence.Type[i] <-"Partial"
  } else{
    metadata_DENV2$Sequence.Type[i] <-"WG"
  }
}

#add updated metadata for vietnam sequences

vietnam_metadata <- read.csv('Data/OUCRU_metadata.csv')

#math to dengue 2

DENV2_vietnam_metadata <- filter(metadata_DENV2, Accession %in% vietnam_metadata$Genome_Accession_number)

minds <- match(DENV2_vietnam_metadata$Accession,vietnam_metadata$Genome_Accession_number)

DENV2_vietnam_metadata$State <- as.matrix(vietnam_metadata$Province..City[minds])

metadata_DENV12<- filter(metadata_DENV2,!Accession %in% DENV2_vietnam_metadata$Accession )

metadata_DENV2$Source <- 'External'
DENV2_vietnam_metadata$Source <- 'Internal'

metadata_DENV2 <- rbind(metadata_DENV2,DENV2_vietnam_metadata)

#omit those with date na

metadata_DENV2 <- metadata_DENV2 %>% 
  filter(!is.na(Date))

#write table for naming process

write.table(metadata_DENV2, file = "Data/metadata_DENV2.txt",row.names=TRUE,sep=",", quote = TRUE)

#plot metadata for sense check 

metadata_DENV2 <- metadata_DENV2 %>% group_by(year = cut(as.Date(Date), "year", start.on.monday = FALSE),Country) 

metadata_DENV2$year <- as.Date(metadata_DENV2$year)

metadata_DENV2_data_plot <- metadata_DENV2 %>%
  group_by(year,Country) %>%
  dplyr :: summarize(count = n())
colnames(metadata_DENV2_data_plot) <- c("year","Country", "count")

ggplot(metadata_DENV2_data_plot, aes(x= as.Date(year), y = count, fill = Country
)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Date of Collection",y="Total Number of Sequences") + theme_bw() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y")

#match to fasta file and name

#make sure to point to the correct directory
path <- "/Users/rhysinward/Documents/Vietnam_nextstrain/Cosmo_manuscript/updated_run/Data/"
fnames<- dir(path)
fnames<- fnames[grep(".fasta",fnames)]
#need to rename file
combinedFname <- paste("sequences.fasta",sep="")
pos <- grep(combinedFname,dir(path),fixed=TRUE)

#read the sequence information file
info <- read.csv(paste(path,"/metadata_DENV2.txt",sep=""))
print(info)
#read file and alter the sequence names 
rootname <-	paste(gsub("\\.fas","",combinedFname),"_beastNames",sep="")
pos <- grep(paste(rootname,".fas",sep=""),dir(),fixed=TRUE)
if (length(pos)==0) {
  
  seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
  taxa <- as.matrix(attributes(seqs)$names)
  taxa <- gsub('.1', '', taxa, fixed = TRUE)
  taxa <- gsub('.2', '', taxa, fixed = TRUE)
  taxa <- gsub(' ', '', taxa, fixed = TRUE)
  genbank_ID <- apply(taxa, 1, getEl, ind=1, sep="\\|")
  minds  <- match(genbank_ID, info$Accession)
  dateTxt <- as.matrix(info$Date[minds])
  sequence_type <- as.matrix(info$Sequence.Type[minds])
  Virus <- 'Dengue_2'
  Strain <- as.matrix(info$Isolate[minds])
  Source <- as.matrix(info$Source[minds])
  decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  country <- as.matrix(info$Country[minds])
  state   <- as.matrix(info$State[minds])
  
  newTaxa <- paste(country,state,genbank_ID,Strain,Virus,sequence_type,Source,dateTxt,decDate,sep="|")
  newTaxa <- gsub(" ","_",newTaxa)
  attributes(seqs)$names <- newTaxa
  write.dna(seqs, file=paste(path,rootname,".fas",sep=""), format="fasta", nbcol=-1, colsep="")
  
  newInfo <- cbind(newTaxa,country,state,Virus,genbank_ID,
                   sequence_type,Source,dateTxt,decDate)
  colnames(newInfo) <- c('Sequence_name',"Country","State","Virus",'genbank_ID',"sequence_type",
                         "Source","Date","decDate")
  write.table(newInfo,file=paste(path,rootname,"_infoTbl.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  write.csv(newInfo,file=paste(path,rootname,"_infoTbl.csv",sep=""),row.names = FALSE)
} else {
  print("Already done renaming")
}

#remove sequences with no associated metadata (removed from earlier process)

seq <- read.fasta("Data/sequencesta_beastNames.fas")
seq_name <- as.data.frame(as.matrix(attributes(seq)$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name$V1),'|',fixed = TRUE)))

#colnames(taxa_split_spatial_temporal) <- c('country','state','strain','virus_id','virus','gene_type',"source",'date','decimal_date')

remove <- filter(seq_name,grepl('NA|NA|NA',V1,fixed = TRUE))

species.to.remove <- remove$V1

vec.names<-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which(! vec.names %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/dengue_denv2.fasta")

