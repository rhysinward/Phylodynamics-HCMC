#Plot case counts
#How do i want to do this? 
#Could try out a few different plots

#Want to have case counts + serotype information with TMRCAs below

#Plot with TMRCA of detected transmission lineages + case counts coloured to indicate denv1/2 +
#panel b with proportions

#Clean Console 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#Load packages 
library(tidyverse)
library(dplyr)
library(lubridate) # for working with dates
library(ggplot2)  # for creating graphs
library(scales)   # to access breaks/formatting functions
library(gridExtra) # for arranging plots
library(data.table)
library(ggforestplot)
library(forcats)

setwd("/Users/rhysinward/Documents/Vietnam_nextstrain/Cosmo_manuscript/updated_run")

#load data

cases <- read.csv('data/denv_case_counts.csv')

# covert date to Date class
cases$YEAR <- as.Date(ISOdate(cases$YEAR, 1, 1))

#convert to long format
long <- melt(setDT(cases), id.vars = c('YEAR'), variable.name = "DENV")

#get proportion of each 

proportion <- long %>%
  group_by(YEAR) %>%
  filter(DENV != 'Total.of.DENV.positive.cases') %>%
  mutate(proportion = value/sum(value)) 

#plot ggplot

#geombar of total cases + line of proportion of cases

#need to get case counts from Tuyen 

long %>%
  filter(DENV == 'Total.of.DENV.positive.cases') %>%
  ggplot( aes(x=YEAR, y=value, group=DENV, color=DENV)) +
  geom_line() + theme_bw() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

#This is several options for plot think I would rather have one plot with cases, one with proportion and one with TMRCA

#geom line for case number
#area for proportion
#tmrca - timeline
coeff <- 80.94

ggplot(proportion, aes(x = YEAR, y = proportion, fill = DENV)) +
  #geom_area(position = "fill", colour = "black", linewidth = .2, alpha = .4) +
  geom_col(aes(y = value/coeff, x = YEAR, fill = DENV),position="dodge", alpha = 1) +
  theme_bw() + labs(x = "Year", y = 'Proportion') +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    # Features of the first axis
    name = "First Axis",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*coeff, name="Second Axis")) 

#area
#get uniform colour scheme for all plots
#remove grid stucture of plot and add gray vs light coloums 
levels(proportion$DENV)
#proportion$DENV= 

proportion$DENV <- factor(proportion$DENV, levels = c('DENV4','DENV3','DENV1','DENV2'))

area <- ggplot(proportion, aes(x = YEAR, y = proportion, fill = DENV)) +
  geom_area(position = "fill", colour = "black", linewidth = .2) +
  theme_bw() + labs(x = "Year", y = 'Proportion of DENV Serotype') +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_fill_manual(name='DENV Serotype',
                    values = c('#C2AFF0', '#A6B07E','#BCA371','#C97064'),
                    breaks=c('DENV1','DENV2','DENV3','DENV4'),
                    labels=c('DENV-1','DENV-2','DENV-3','DENV-4')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",),
        text = element_text(size = 40))

#cases 

cases <- long %>%
  filter(DENV == 'Total.of.DENV.positive.cases') %>%
  ggplot(aes(x=YEAR, y=value, group=DENV, color=DENV)) +
  geom_line(colour = 'black') +
  theme_bw() + labs(x = "Year", y = 'Number of Confirmed DENV Cases') +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40))

#tmrca 

#want this on same scale 2008 - 2022 - want to have a dot or something for central value 
#with the CI presndted by a dotted line? remove horizontal on plot

### Load packages
library("vistime")
library("tidyverse")
library("RColorBrewer")
library("scales")
library("cowplot")
library('patchwork')

### Generate input data
tmrca_clade <- data.frame(clade = c("Clade A", "Clade B", "Clade C"),
                         start = c("2015-06-01", "2018-07-01", "2016-09-01"),
                         end = c("2017-03-01", "2021-01-01", "2018-08-01"),
                         tmrca = c("2016-04-01", "2019-12-01", "2017-08-01"))

### Preparing for plot
## Colours
# Define number of colours needed
cols_n <- as.numeric(n_distinct(tmrca_clade$clade))

# Check we have enough colours
ifelse(cols_n>12, 
       "More than 12 clades - not enough colours!",
       "12 clades or fewer - using Set3 colours")

# Select the colours from Set3 palette in RColorBrewer
cols_to_use <- brewer.pal(n = cols_n, name = "Set3")

# Create mapping of colours to clade
col_clade_mapping <- data.frame(clade=unique(c(as.character(tmrca_clade$clade))), color=cols_to_use)

# merge in the mapping to the df
clade_moves_2 <- merge(tmrca_clade,
                      col_clade_mapping,
                      by="clade",
                      all.x=T,all.y=T) %>%
  select(clade, tmrca, start, end, color) 
clade_moves_2

## Extract clade dates
clade_dates <- clade_moves_2 %>%
  select(clade, tmrca,start, end, color) %>%
  distinct(tmrca, .keep_all=TRUE)


### Plotting
# Produce the basic plot
plot_data <- gg_vistime(data = clade_dates,
                        col.group = "clade", # Each row will be a patient
                        col.event = "clade", # Rows will be coloured by the ward
                        show_labels = FALSE, # Remove labels indicating the ward
                        linewidth = 10)

# Tweak the plot

# define the alternating colors for the background
colors <- c("#F0F0F0", "#E0E0E0")

# plot the data with a background color for each year

plot_data <- plot_data + theme_bw() +
  ggplot2::theme(
    plot.title = element_text(size=14),
    axis.text.x = element_text(size = 15, color = "black", angle = 30, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15, color = "black")) +
  scale_x_datetime(breaks = breaks_width("1 year"), labels = date_format("%Y")) +
  theme_minimal() +
  scale_color_manual(values = c('#64A7C9','#7E74B0','#71BCA3')) +
  theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40)) +
  labs(x = "Year")

plot_data

# Adding tmrca date
plot_data <- plot_data +
  annotate("point", x = as.POSIXct(tmrca_clade[1,4]), y = 5, size = 8, colour = "black") +
  annotate("point", x = as.POSIXct(tmrca_clade[2,4]), y = 3, size = 8, colour = "black") +
  annotate("point", x = as.POSIXct(tmrca_clade[3,4]), y = 1, size = 8, colour = "black")
plot_data

#arrange plots

arranged_plot <- (cases | area) / plot_data + plot_annotation(tag_levels ="A")

### Save plots
ggplot2::ggsave(arranged_plot, file = "arranged.png", height=20, width=35)



