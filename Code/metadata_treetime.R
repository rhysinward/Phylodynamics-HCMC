#Metadata for Treetime

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

library(tidyverse)
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
library(tidygeocoder)
library(treeio)


#change sequence name to match data - iqtree2 has changed () --> _

treemmer_seq <- seqinr :: read.fasta("Results/denv2_E_gene_cosmo_cleaned.fasta")
treemmer_seq_name <- as.data.frame(as.matrix(attributes(treemmer_seq)$names))
treemmer_seq_name$V1 <- gsub("(", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub(")", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub(",", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub("*", "_", treemmer_seq_name$V1, fixed = TRUE)

write.fasta(sequences=treemmer_seq, names=treemmer_seq_name$V1, file.out="Results/denv2_E_gene_cosmo_cleaned.fasta")

#denv_2

treemmer_seq <- seqinr :: read.fasta("Results/denv2_E_gene_cosmo_cleaned.fasta")
treemmer_seq_name <- as.data.frame(as.matrix(attributes(treemmer_seq)$names))
taxa_split_spatial_temporal <- data.frame(do.call('rbind',strsplit(as.character(treemmer_seq_name$V1),'|',fixed = TRUE)))

taxa_split_spatial_temporal$X1[taxa_split_spatial_temporal$X1 == "Viet_Nam"] <- "Vietnam"

test <- select(taxa_split_spatial_temporal, c(1,2,5,8,9))

at_longs <- test %>%
  geocode(X1, method = 'osm', lat = latitude , long = longitude)

taxa_split_spatial_temporal <- select(taxa_split_spatial_temporal, c(1,2,6,8,9))

taxa_split_spatial_temporal$latitude <- at_longs$latitude
taxa_split_spatial_temporal$longitude <- at_longs$longitude

taxa_split_spatial_temporal$region <- 'SEA'

taxa_split_spatial_temporal$X8 <- as.Date(taxa_split_spatial_temporal$X8)

taxa_split_spatial_temporal$num_date <- year(taxa_split_spatial_temporal$X8)

taxa_split_spatial_temporal <- cbind(taxa_split_spatial_temporal,treemmer_seq_name)

taxa_split_spatial_temporal <- taxa_split_spatial_temporal %>% 
  dplyr :: rename(
    country = X1,
    state = X2,
    genome_type = X6,
    cal_date = X8,
    date = X9,
    strain = V1
  )

taxa_split_spatial_temporal <- select(taxa_split_spatial_temporal, c(10,1,2,3,4,5,6,7,8,9))

write.csv(taxa_split_spatial_temporal,file="Results/metadata_denv2_cosmo.csv",row.names = FALSE,quote=FALSE)


#same but for DENV-2 all serotypes 

treemmer_seq <- seqinr :: read.fasta("Results/denv2_E_gene_cleaned.fasta")
treemmer_seq_name <- as.data.frame(as.matrix(attributes(treemmer_seq)$names))
treemmer_seq_name$V1 <- gsub("(", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub(")", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub(",", "_", treemmer_seq_name$V1, fixed = TRUE)
treemmer_seq_name$V1 <- gsub("*", "_", treemmer_seq_name$V1, fixed = TRUE)

write.fasta(sequences=treemmer_seq, names=treemmer_seq_name$V1, file.out="Results/denv2_E_gene_cleaned.fasta")

#denv_2

treemmer_seq <- seqinr :: read.fasta("Results/denv2_E_gene_cleaned.fasta")
treemmer_seq_name <- as.data.frame(as.matrix(attributes(treemmer_seq)$names))
taxa_split_spatial_temporal <- data.frame(do.call('rbind',strsplit(as.character(treemmer_seq_name$V1),'|',fixed = TRUE)))

taxa_split_spatial_temporal$X1[taxa_split_spatial_temporal$X1 == "Viet_Nam"] <- "Vietnam"

test <- select(taxa_split_spatial_temporal, c(1,2,5,8,9))

at_longs <- test %>%
  geocode(X1, method = 'osm', lat = latitude , long = longitude)

taxa_split_spatial_temporal <- select(taxa_split_spatial_temporal, c(1,2,6,8,9))

taxa_split_spatial_temporal$latitude <- at_longs$latitude
taxa_split_spatial_temporal$longitude <- at_longs$longitude

taxa_split_spatial_temporal$region <- 'SEA'

taxa_split_spatial_temporal$X8 <- as.Date(taxa_split_spatial_temporal$X8)

taxa_split_spatial_temporal$num_date <- year(taxa_split_spatial_temporal$X8)

taxa_split_spatial_temporal <- cbind(taxa_split_spatial_temporal,treemmer_seq_name)

taxa_split_spatial_temporal <- taxa_split_spatial_temporal %>% 
  dplyr :: rename(
    country = X1,
    state = X2,
    genome_type = X6,
    cal_date = X8,
    date = X9,
    strain = V1
  )

taxa_split_spatial_temporal <- select(taxa_split_spatial_temporal, c(10,1,2,3,4,5,6,7,8,9))

write.csv(taxa_split_spatial_temporal,file="Results/metadata_denv2_all.csv",row.names = FALSE,quote=FALSE)


