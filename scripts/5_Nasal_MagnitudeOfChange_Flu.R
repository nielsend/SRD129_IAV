#######################################################################
#SRD129 16S - Magnitude of Change in Nasal Microbiota


#Purpose: Plot the F-statistic from PERMANOVA pairwise comparisons of Control and Flu groups, over time (displays the magnitude of change in the nasal bacterial community structure of the Flu group relative to control)

#Files needed:
#Flu_MagnitudeOfChange.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library(ggrepel)


#######################################################################

#The F-values were obtained from 'Flu.pairwisecomparisons.csv' data generated in the "Beta diversity" section
#Flu_MagnitudeOfChange.csv was created from Flu.pairwisecomparisons.csv
#To make the Flu_MagnitudeOfChange.csv file, open Flu.pairwisecomparisons.csv file in excel, 
#copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. Add "Day" and "Treatment" columns and save as "Flu_MagnitudeOfChange.csv".
View(Flu.pairwisecomparisons)

flu.mag <- read.csv("./Flu_MagnitudeOfChange.csv")
colnames(flu.mag) <- c("F.Model","R2","p.value","p.adjusted","p.adjusted2","Day","Treatment")
class(flu.mag)
flu.mag$Day <- factor(flu.mag$Day) #Encode "Day" as a factor
flu.mag$Day = factor(flu.mag$Day, levels=c("D0", "D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))

flu.mag2 <- ggplot(data=flu.mag, aes(x=Day, y=F.Model, group=Treatment)) +
  geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("#00BFC4")) +
  scale_color_manual(values=c("#00BFC4")) +
  geom_label_repel(aes(label=p.adjusted2), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50') +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group")
flu.mag2

#Save 'flu.mag2' as a .tiff for publication, 500dpi
ggsave("Figure_3.tiff", plot=flu.mag2, width = 10, height = 5, dpi = 500, units =c("in"))
