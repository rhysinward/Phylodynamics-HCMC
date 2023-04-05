#Create trees for cosmo manuscript 
#This will create 3 seperate trees 
#1 ML with all of DENV-2 with genotype labeled
#1 Cosmo ML with clades with OURCU Vietnam Sequences zoomed into 
#1 time-scaled cosmo tree with TMRCA for these clades 

###want to add black boarder to the tips + increase branch width 
###want to have a plot showing what Vi plot shows + add tmrca as a second panal

#Clean Console 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#Load packages 
library(phangorn)
library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)
library(plotly)
library(seqinr)
library(dispRity)
library(rstudioapi)
library(fBasics)
library(grDevices)

#SWD

setwd("/Users/rhysinward/Documents/Vietnam_nextstrain/Cosmo_manuscript/updated_run")

#Load in cosmo tree

tree <- read.tree('Results/denv2_E_gene_cosmo_cleaned.fasta.treefile')

#Get metadata from tip name

metadata_df <- as.data.frame(tree[["tip.label"]])

metadata <- data.frame(do.call('rbind',strsplit(as.character(metadata_df$`tree[["tip.label"]]`),'|',fixed = TRUE)))

metadata <- cbind(metadata_df,metadata)

metadata <- select(metadata, c(1,2,4,7,8,9,10))

metadata <- metadata %>% 
  dplyr :: rename(
    country = X1,
    accession = X3,
    genome_type = X6,
    date = X8,
    dec_date = X9,
    strain = `tree[["tip.label"]]`
  )

#metadata[1,6] <- '2006-06-15'
#metadata[1,2] <- 'Viet_Nam'
#metadata[2,6] <- '2006-06-15'
#metadata[2,2] <- 'Viet_Nam'

write.csv(metadata,file="metadata_denv2_cosmo.csv",row.names = FALSE)

#manually added info of added sequences to this csv

#metadata processing done

x <- read.csv("metadata_denv2_cosmo.csv")

#want to edit change the country information to bordering, other, Vietnam, vietnam_sequenced

#remove ' from tree name

x$country <-gsub("'","",as.character(x$country))

#create lookup table for OURCU sequences

look_up <- c("OP090394","OP090395", "OP090396","OP090397","OP090390","OP090391","OP090392",
             "OP090388","OP090393","OP090389","OP458511","OP458514","OP458516","OP458515","OP458512","OP458517",
             "OP458513","OQ028208","OQ028206","OQ028207","OQ028212","OQ028211","OQ028210","OQ028209","OQ028205","OQ028213",
             "OQ028229","OQ028230","OQ028231","OQ028215","OQ028216","OQ028217","OQ028218","OQ028219","OQ028220","OQ028221",
             "OQ028222","OQ028223","OQ028224","OQ028225","OQ028226","OQ028227","OQ028228","OQ028232","OQ028214")

#add location info

taxa_1 <- x
taxa_1$Location <- NA

for (i in 1:nrow(taxa_1)){
  if (taxa_1$country[i] == 'China'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Laos'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Cambodia'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Viet_Nam'){
    taxa_1$Location[i] <- 'Other Vietnam'
  } else {
    taxa_1$Location[i] <- 'Other Southern and SEA'
  }
}

taxa_2 <- filter(taxa_1,accession %in% look_up)
taxa_2$Location = ('OUCRU-HCMC Vietnam')
taxa_3 <- filter(taxa_1,!accession %in% look_up)

taxa_4 <- rbind(taxa_2,taxa_3)

#plot tree
p <- ggtree(tree,size=1) %<+% taxa_4

p1 <- p +
  geom_tippoint(aes(color=Location), size=6, alpha=.75) +
  scale_color_manual(values = c('#C2AFF0','#BCA371','#A6B07E','#C97064')) +
  theme(legend.position="right") +
  geom_text2(aes(subset = !isTip, label=label), # subset=!isTip
             size = 3.5,
             color = "black",
             hjust = 1, 
             vjust = -1.5
  )

p1 <- p1 + geom_text2(aes(subset = !isTip, label=label))

#p1 <- p1 +   geom_strip('Cambodia|NA|OL414742|100-0153_2019-06-23|Dengue_2|WG|2019-06-15|2019.45205479452',
#                  'Cambodia|NA|OL414747|100-0034_2019-08-20|Dengue_2|WG|2019-08-15|2019.61917808219',
#                  barsize=3, color='red', 
#                  label="associated taxa", offset.text=.1) +
#  geom_strip('Taiwan|NA|MG895160||Dengue_2|E|2016-09-01|2016.66575342466',
#             'Indonesia|Purwokerto|KY709187|PWT-077|Dengue_2|E|2015-06-15|2015.45205479452',
#             barsize=3, color='blue', 
#             label = "another label", offset.text=.1) +
#  geom_strip('Viet_Nam|NA|OP458515|56DX-105|Dengue_2|E|2021-02-22|2021.14246575342',
#             'Indonesia|Yogyakarta|OK180523|MJ1-07-11|Dengue_2|E|2016-07-21|2016.55068493151',
#             barsize=3, color='green', 
#             label = "another label", offset.text=.1)

#p2 <- collapse(p1, node=4720) + 
#  geom_point2(aes(subset=(node==4720)), shape=23, size=5,alpha=.75, fill='black') #example code for collapsing nodes

#View clades with OURCU sequences

metat <- p1$data %>%
  dplyr::inner_join(taxa_4, c('label' = 'strain'))

p3 <- p1 +
  geom_point(data = metat,
             aes(x = x,
                 y = y,
                 colour = Location.y
             ))

Iplot <- plotly::ggplotly(p3) #creates interative tree to get node position

p4 <- p +
  geom_tippoint(aes(color=Location), size=7, alpha=1,shape = 21, colour = 'black' ) +
  geom_tippoint(aes(color=Location), size=7, alpha=.7) +
  scale_color_manual(values = c('#C2AFF0','#A6B07E','#BCA371','#dc5349')) +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 2830)

p5 <- p +
  geom_tippoint(aes(color=Location), size=15, alpha=1,shape = 21, colour = 'black' ) +
  geom_tippoint(aes(color=Location), size=15, alpha=.7) +
  scale_color_manual(values = c('#C2AFF0','#A6B07E','#BCA371','#dc5349')) +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 2830)

p4

#p5 <- collapse(p4, node=3122) + 
#  geom_point2(aes(subset=(node==3122)), shape=23, size=3,alpha=.75, fill='black')

#p6 <- collapse(p5, node=3064) + 
#  geom_point2(aes(subset=(node==3064)), shape=23, size=3,alpha=.75, fill='black')

clade_1 <- viewClade(p5, 4865)
clade_2 <- viewClade(p5, 3769)
clade_3 <- viewClade(p5, 2830)
clade_1_bb <- viewClade(p1, 4865)
clade_2_bb <- viewClade(p1, 3769)
clade_3_bb <- viewClade(p1, 2830)

#Save plot outputs

ggsave("cosmo_ml_tree.png", plot = p4,width = 20, height = 15)
ggsave("cosmo_ml_clade1.png", plot = clade_1,width = 30, height = 10)
ggsave("cosmo_ml_clade2.png", plot = clade_2,width = 30, height = 10)
ggsave("cosmo_ml_clade3.png", plot = clade_3,width = 40, height = 10)
ggsave("cosmo_ml_clade1_bb.png", plot = clade_1_bb,width = 30, height = 10)
ggsave("cosmo_ml_clade2_bb.png", plot = clade_2_bb,width = 30, height = 10)
ggsave("cosmo_ml_clade3_bb.png", plot = clade_3_bb,width = 40, height = 10)

knitr::plot_crop('cosmo_ml_clade1.png')
knitr::plot_crop('cosmo_ml_clade2.png')
knitr::plot_crop('cosmo_ml_clade3.png')
knitr::plot_crop('cosmo_ml_tree.png')


#time-scaled cosmo tree 

#load tree

egtree <- treeio :: read.nexus('Results/timetree_cosmo.nexus')

#egtree$edge.length <- egtree$edge.length/0.00086 #need to divide edge length by mutational rate

tt <- ggtree(rescale_tree(egtree, branch.length = 'clock'), mrsd = max(x$date),size=1) +
  coord_cartesian(clip = 'off') +  theme_tree2() 

tt <- tt %<+% taxa_4

age <- tree.age(egtree, order = "past", fossil = FALSE) #this gives the age of all 

age$ages <- max(x$dec_date)-age$ages

tt_plot <- tt +
  geom_tippoint(aes(color=Location), size=5, alpha=1,shape = 21, colour = 'black' ) +
  geom_tippoint(aes(color=Location), size=5, alpha=.7) +
  scale_color_manual(values = c('#C2AFF0','#BCA371','#A6B07E','#dc5349')) +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 2820) +
  theme(text = element_text(size = 40)) +
  geom_vline(xintercept=c(1990,2000,2010,2020), linetype="dotted") 

#get node for the OURCU sequences

tt_plot <- tt +
  geom_tippoint(aes(color=Location), size=6, alpha=.75) +
  scale_color_brewer("Location", palette="Spectral") +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 1550) +
  theme(text = element_text(size = 40)) +
  geom_vline(xintercept=c(1990,2000,2010,2020), linetype="dotted") 

tt_plot <- tt +
  geom_tippoint(aes(color=Location), size=6, alpha=.75) +
  scale_color_brewer("Location", palette="Spectral") +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 3000) +
  theme(text = element_text(size = 40)) +
  geom_vline(xintercept=c(1990,2000,2010,2020), linetype="dotted") +
  geom_vline(xintercept=c(2016.689,2020.176,2017.563), linetype="dashed") 

metat_1 <- tt_plot$data %>%
  dplyr::inner_join(taxa_4, c('label' = 'strain'))

tt_plot_1 <- tt_plot +
  geom_point(data = metat_1,
             aes(x = x,
                 y = y,
                 colour = label
             ))

plotly::ggplotly(tt_plot_1) #creates interative tree to get node position

#get node names from icytree as ggtree changes the names
viewClade(tt_plot_1, 1813)  
viewClade(tt_plot_1, 2079) 
viewClade(tt_plot_1, 1947) 

ggsave("Results/time_scaled_comso.png", plot = tt_plot,width = 15, height = 10)

#now we want to process file with all sequences

#Match the treetime names to the fasta file to get only the sequences we are interested
#and run these through tree-time

#Load in cosmo tree

tree <- read.tree('Results/denv2_E_gene_cleaned.fasta.treefile')

#Get metadata from tip name

metadata_df <- as.data.frame(tree[["tip.label"]])

metadata <- data.frame(do.call('rbind',strsplit(as.character(metadata_df$`tree[["tip.label"]]`),'|',fixed = TRUE)))

metadata <- cbind(metadata_df,metadata)

metadata <- select(metadata, c(1,2,4,7,8,9,10))

metadata <- metadata %>% 
  dplyr :: rename(
    country = X1,
    accession = X3,
    genome_type = X6,
    date = X8,
    dec_date = X9,
    strain = `tree[["tip.label"]]`
  )

#metadata[1,6] <- '2006-06-15'
#metadata[1,2] <- 'Viet_Nam'
#metadata[2,6] <- '2006-06-15'
#metadata[2,2] <- 'Viet_Nam'

write.csv(metadata,file="metadata_denv2_all.csv",row.names = FALSE)

#manually added info of added sequences to this csv

#metadata processing done

x <- read.csv("metadata_denv2_all.csv")

#want to edit change the country information to bordering, other, Vietnam, vietnam_sequenced

#remove ' from tree name

x$country <-gsub("'","",as.character(x$country))

#create lookup table for OURCU sequences

look_up <- c("OP090394","OP090395", "OP090396","OP090397","OP090390","OP090391","OP090392",
             "OP090388","OP090393","OP090389","OP458511","OP458514","OP458516","OP458515","OP458512","OP458517",
             "OP458513","OQ028208","OQ028206","OQ028207","OQ028212","OQ028211","OQ028210","OQ028209","OQ028205","OQ028213",
             "OQ028229","OQ028230","OQ028231","OQ028215","OQ028216","OQ028217","OQ028218","OQ028219","OQ028220","OQ028221",
             "OQ028222","OQ028223","OQ028224","OQ028225","OQ028226","OQ028227","OQ028228","OQ028232","OQ028214")

#add location info

taxa_1 <- x
taxa_1$Location <- NA

for (i in 1:nrow(taxa_1)){
  if (taxa_1$country[i] == 'China'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Laos'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Cambodia'){
    taxa_1$Location[i] <- 'Border'
  } else if (taxa_1$country[i] == 'Viet_Nam'){
    taxa_1$Location[i] <- 'Other Vietnam'
  } else {
    taxa_1$Location[i] <- 'Other Southern and SEA'
  }
}

taxa_2 <- filter(taxa_1,accession %in% look_up)
taxa_2$Location = ('OUCRU-HCMC Vietnam')
taxa_3 <- filter(taxa_1,!accession %in% look_up)

taxa_4 <- rbind(taxa_2,taxa_3)

#load genome detective output and merge with metadata

serotype <- read.csv('Results/metadata_genome_detective.csv')
seq_name <- data.frame(do.call('rbind',strsplit(as.character(serotype$name),'|',fixed = TRUE)))
serotype$accession <- seq_name$X3

serotype <- select(serotype,-c(1))

taxa_4 <- taxa_4 %>% left_join( serotype, 
                                by=c('accession'='accession'))

#plot tree

#plot tree

p <- ggtree(tree, size = 1) %<+% taxa_4

p1 <- p +
  geom_tippoint(aes(color=Location), size=6, alpha=.75) +
  theme(legend.position="none") 

p4 <- p +
  geom_tippoint(aes(color=Location), size=5, alpha=1,shape = 21, colour = 'black' ) +
  geom_tippoint(aes(color=Location), size=5, alpha=.7) +
  scale_color_manual(values = c('#C2AFF0','#A6B07E','#BCA371','#dc5349')) +
  theme(legend.position="none") +
  ggplot2::ylim(-0.5, 4500)

p1 <- collapse(p1, node=2931) + 
  geom_point2(aes(subset=(node==2931)), shape=23, size=3,alpha=.75, fill='black') #example code for collapsing nodes
p1 <- collapse(p1, node=2902) + 
  geom_point2(aes(subset=(node==2902)), shape=23, size=3,alpha=.75, fill='black') #example code for collapsing nodes



#View clades get names for plot sequences

metat <- p1$data %>%
  dplyr::inner_join(taxa_4, c('label' = 'strain'))

p3 <- p1 +
  geom_point(data = metat,
             aes(x = x,
                 y = y,
                 colour = label))

#plotly::ggplotly(p3) #creates interative tree to get node position

#plot editiing 
p3 <- p4 + geom_strip('China|NA|MK543450|HNQY20180091|Dengue_2|WG|External|2018-09-27|2018.73698630137',
                      'India|NA|MZ277527|14441_Jharkhand_India_2018|Dengue_2|E|External|2018-04-25|2018.31232876712',
                      barsize=3, color='#36454F', 
                      offset.text=.1) +
  geom_strip('China|NA|KX577643|918-06|Dengue_2|E|External|2015-11-01|2015.83287671233',
             'Singapore|NA|JN030344|SG_EHI_D2/31961Y10|Dengue_2|E|External|2010-06-15|2010.45205479452',
             barsize=3, color='#964B00', 
             offset.text=.1) +
  geom_strip('China|Guangzhou|MN913509|817/2018|Dengue_2|E|External|2018-06-15|2018.45205479452',
             'India|NA|MH594922|NIV_1736862_Nashik_India_2017|Dengue_2|E|External|2017-06-15|2017.45205479452',
             barsize=3, color='#FADA5E', 
             offset.text=.1) +
  geom_strip('Myanmar|NA|KX357976|MandalayA12|Dengue_2|E|External|2015-06-15|2015.45205479452',
             'India|NA|OP809583|THSTI-TRC-DENV2-12|Dengue_2|WG|External|2019-06-15|2019.45205479452',
             barsize=3, color='#AFE1AF',
             offset.text=.1)

ggsave("cosmo_tt_tree.png", plot = p4,width = 15, height = 10)


#p1 <- collapse(p1, node=2931) + 
#    geom_point2(aes(subset=(node==2931)), shape=23, size=3,alpha=.75, fill='black') #example code for collapsing nodes
#p1 <- collapse(p1, node=2902) + 
#  geom_point2(aes(subset=(node==2902)), shape=23, size=3,alpha=.75, fill='black') #example code for collapsing nodes

#metat <- p1$data %>%
#  dplyr::inner_join(taxa_4, c('label' = 'name'))

#p3 <- p1 +
#  geom_point(data = metat,
#             aes(x = x,
#                 y = y,
#                 colour = label))

#plotly::ggplotly(p3) #creates interative tree to get node position

#Save plot outputs

#ggsave("cosmo_tt_tree.png", plot = p1,width = 15, height = 10)

