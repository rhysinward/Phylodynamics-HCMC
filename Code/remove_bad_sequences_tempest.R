#Remove bad sequences through TempEst

###Code Start###

setwd("C:/Users/rhysi/OneDrive/Documents/Oxford/trail_run_for_github")

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

tempest <- read.csv('results/DENV2_EG_tempest.txt',sep='')

DENV2_EG_filtered <- filter(tempest, residual <= sd(tempest$residual) & residual >=
                              -sd(tempest$residual))

DENV2 <- read.fasta("results/aligned.fasta")

DENV2_EG <- DENV2[names(DENV2) %in% DENV2_EG_filtered$tip]

write.dna(DENV2_EG, file=("results/denv2_E_gene_QC.fasta"), format="fasta",
          nbcol=-1, colsep="")

#DENV 2 WG

tempest <- read.csv('Results/DENV2_WG_tempest.txt',sep='')

DENV2_WG_filtered <- filter(tempest, residual <= sd(tempest$residual) & residual >=
                              -sd(tempest$residual))

DENV2 <- read.fasta("Results/denv2_WG.fasta")

DENV2_WG <- DENV2[names(DENV2) %in% DENV2_WG_filtered$tip]

write.dna(DENV2_WG, file=("results/denv2_WG_QC.fasta"), format="fasta",
          nbcol=-1, colsep="")
