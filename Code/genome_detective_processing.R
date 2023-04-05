#create data for genome detective and process the outputs 

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

#Load sequences

seq <- read.fasta("Data/dengue_denv2.fasta")
seq_name <- as.data.frame(as.matrix(attributes(seq)$names))

#colnames(taxa_split_spatial_temporal) <- c('country','state','strain','virus_id','virus','gene_type',"source",'date','decimal_date')
#Split into fasta files of 2000 sequences max

keep <- seq_name[1:2000,]

vec.names<-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which( vec.names %in%  keep)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/dengue_denv2_part1.fasta")

keep <- seq_name[2001:4000,]

vec.names<-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which( vec.names %in%  keep)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/dengue_denv2_part2.fasta")

keep <- seq_name[4001:5836,]

vec.names<-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which( vec.names %in%  keep)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/dengue_denv2_part3.fasta")

#Process outputs of Genome Detective

part_1 <- read.csv('Data/results.csv')
part_2 <- read.csv('Data/results_(1).csv')
part_3 <- read.csv('Data/results_(2).csv')

serotype_2 <- rbind(part_1,part_2,part_3)

seq <- read.fasta("Data/dengue_denv2.fasta")
seq_name <- as.data.frame(as.matrix(attributes(seq)$names))
serotype_2$name <- seq_name$V1

vec.tokeep <-which(serotype_2$species ==  'Dengue virus type 2')

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/dengue_denv2_unaligned.fasta")
write.csv(serotype_2,file="results/metadata_genome_detective.csv",row.names = FALSE)

