#Split into WG and E gene and remove sequences with a high number of N's

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
library(Biostrings)

#get E gene # this maybe redundant

full_alignment <- read.dna("Results/aligned.fasta", format="fasta", as.matrix=TRUE)
gene_wanted <- full_alignment[,937:2421]
write.dna(gene_wanted, file="Results/denv2_E_gene.fasta", format="fasta", nbcol=-1, colsep="")

#remove sequences with high proportion of N

E_gene <- read.fasta("Results/denv2_E_gene.fasta")

seqs <- rapply(E_gene,function(x) ifelse(0.05*sum(table(x)) > sum(x == 'n'),'good',x), how = "unlist")

vec.tokeep <-which(seqs ==  'good')

write.fasta(sequences=E_gene[vec.tokeep], names=names(E_gene)[vec.tokeep], file.out="Results/denv2_E_gene.fasta")

#get cosmo

seq <- read.fasta("Results/denv2_E_gene.fasta")
seq_name <- as.data.frame(as.matrix(attributes(seq)$names))
gd <- read.csv('Results/metadata_genome_detective.csv')

remove <- filter(gd, type != 'DENV-2 Genotype II - Cosmopolitan')

species.to.remove <- remove$name

vec.names <-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which(! vec.names %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Results/denv2_E_gene_cosmo.fasta")
