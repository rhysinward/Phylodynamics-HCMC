#Remove bad sequences through TempEst

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
library(ShortRead) #BiocManager::install("ShortRead")

#DENV 2 EG

tempest <- read.csv('Results/DENV2_tempest.txt',sep='')
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(tempest$tip),'|',fixed = TRUE)))
tempest$accession <- taxa_split$X3

DENV2_EG_filtered <- tempest %>%
  mutate(year = floor(date)) %>%
  group_by(year) %>%
  filter(distance <= (sd(tempest$distance)*4 + mean(tempest$distance))
         & distance >= (mean(tempest$distance) -sd(tempest$distance)*4))

DENV2 <- read.fasta("Results/denv2_E_gene.fasta")
seq_name <- as.data.frame(as.matrix(attributes(DENV2)$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name$V1),'|',fixed = TRUE)))

vec.tokeep <-which(taxa_split$X3 %in%  DENV2_EG_filtered$accession)

length(vec.tokeep)

write.fasta(sequences=DENV2[vec.tokeep], names=names(DENV2)[vec.tokeep], file.out="Results/denv2_E_gene_cleaned.fasta")


#DENV 2 Cosmo

tempest <- read.csv('Results/DENV2_cosmo_tempest.txt',sep='')
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(tempest$tip),'|',fixed = TRUE)))
tempest$accession <- taxa_split$X3

DENV2_EG_filtered <- tempest %>%
  mutate(year = floor(date)) %>%
  group_by(year) %>%
  filter(distance <= (sd(tempest$distance)*4 + mean(tempest$distance))
         & distance >= (mean(tempest$distance) -sd(tempest$distance)*4))

DENV2 <- read.fasta("Results/denv2_E_gene_cosmo.fasta")
seq_name <- as.data.frame(as.matrix(attributes(DENV2)$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name$V1),'|',fixed = TRUE)))

vec.tokeep <-which(taxa_split$X3 %in%  DENV2_EG_filtered$accession)

length(vec.tokeep)

write.fasta(sequences=DENV2[vec.tokeep], names=names(DENV2)[vec.tokeep], file.out="Results/denv2_E_gene_cosmo_cleaned.fasta")


