#######################################################################
#SRD129 Nasal Magnitude of Change
#Kathy Mou

#NOTES: 
#This code plots the magnitude of change for nasal microbiota beta diversity

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on desktop or network drive)
#Mac
setwd("~/Desktop/SRD129/SRD129_2000singletons")

#Load library packages
library(vegan)
library(ggplot2)
#ggplot:  fill = color in boxes
#         color = colored outline
#aes: aesthetics declared
#alt - to get <- symbol
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
install.packages("ggrepel")
library(ggrepel)

#Load image file
load("SRD129_NasalTonsil_MagnitudeOfChange2.RData")





#Save image file
save.image(file="SRD129_NasalTonsil_MagnitudeOfChange2.RData")

######################################################################

#The F-values were obtained from SRD129_Nasal_SelectPairwiseAdonisFluControlNoDNEG.csv generated in the "Beta diversity" section
#SRD129_Nasal_FluMagnitudeofChange.csv was created from SRD129_Nasal_SelectPairwiseAdonisFluControlNoDNEG.csv data that was rearranged
#To make the SRD129_Nasal_FluMagnitudeofChange.csv file, open SRD129_Nasal_SelectPairwiseAdonisFluControlNoDNEG.csv file 
#(created from SRD129_alpha_beta_diversity_Flu.R) in excel, copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. 
#Add "Day" and "Treatment" columns and save as "SRD129_Nasal_FluMagnitudeofChange.csv".

#Nasal
nasalmag <- read.csv("SRD129_Nasal_FluMagnitudeofChange.csv",  header = TRUE, sep = "")
class(nasalmag) #data.frame
nasalmag$Day <- factor(nasalmag$Day)
nasalmag
nasal <- ggplot(data=nasalmag, aes(x=Day, y=F.Model, group=Treatment)) +
  geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("#00BFC4")) +
  scale_color_manual(values=c("#00BFC4")) +  
  geom_label_repel(aes(label=p.adjusted2), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50') +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size=14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14)) +
  labs(color="Treatment group")
nasal

#Save 'nasal' as a .tiff for publication, 500dpi
ggsave("Figure_3.tiff", plot=nasal, width = 10, height = 5, dpi = 500, units =c("in"))
