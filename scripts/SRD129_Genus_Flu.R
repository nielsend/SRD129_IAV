############################################################
#SRD129 IAV, Control total genera
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs
#Only analyzed data without bad samples

#Set working directory
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("~/Desktop/SRD129/SRD129_2000singletons")
setwd("C:/Users/Kathy.Mou/Desktop/SRD129/SRD129_2000singletons")

#Clear workspace and load necessary packages
rm(list=ls())

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")
library(data.table)

#Load saved image
load("SRD129.Genus2.RData")



#Save image
save.image(file="SRD129.Genus2.RData")

############################################################

##No bad samples metadata

otu <- import_mothur(mothur_shared_file = 'SRD129.outsingletons.abund.opti_mcc.0.03.subsample.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'SRD129.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('SRD129.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('SRD129metadata06292018nobad.csv', header = TRUE, sep = ",")
head(meta)

#############NASAL#############
#Subset each "All" group from fobar2.gather and save as individual csv files.
#For example:
D0_Control_NW.genus <- subset(fobar2.gather, All=="D0_Control_NW") #EDIT THIS
write.csv(D0_Control_NW.genus, file="D0_Control_NW.genus.csv")  #EDIT THIS

#Calculate the total % percent abundance of each genera on each sample (I used JMP to do this) 
#and save results in a spreadsheet editor such as Excel (see D0_Control_NW.genus.xlsx for an example)
#Since we are only interested in genera that are above 2% abundance, 
#calculate total percentage of all other genera that are less than 2% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "Nasal genus.csv" or "Tonsil genus.csv". 
#Create the following columns: Day, Treatment group, Tissue, Percent Abundance, and Genus. 
#Copy the list of genera and their percent abundances from each of the individual Excel files to the respective "Nasal genus.csv" or "Tonsil genus.csv" spreadsheet.
#Fill in the other columns manually (Day, Treatment Group, Tissue). 
#You should have a file similar to SRD129_Nasal_GenusPercentAbundance.csv. Continue to the next step.
nasalgen <- read.table('SRD129_Nasal_GenusPercentAbundance.csv', header = TRUE, sep = ",")
head(nasalgen)

#Nasalgen.2 (flu and control only, more than 1% genera, remove DNEG12, DNEG6)
#nasalgen.2<- nasalgen[!grepl("PRRSV", nasalgen$Treatment),]
#nasalgen.2<- nasalgen.2[!grepl("BB", nasalgen.2$Treatment),]
nasalgen.2 <- nasalgen.2[!grepl("DNEG12", nasalgen.2$Day),]
nasalgen.2 <- nasalgen.2[!grepl("DNEG6", nasalgen.2$Day),]
unique(nasalgen.2$Treatment) #Control IAV
unique(nasalgen.2$Day) #D0  D1  D3  D7  D10 D14 D21 D28 D36 D42
nasalgen.2$Day = factor(nasalgen.2$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
nasalgen.2$More.than.2=as.character(nasalgen.2$Genus)
str(nasalgen.2$More.than.2) #Compactly display the internal structure of an R object,
nasalgen.2$More.than.2[nasalgen.2$Percent.abundance<1]<-"Other"
write.csv(nasalgen.2, file = "SRD129_Nasal_FluControlNoDNEGGenusPercentAbundanceAbove1percent.csv")

#To make sure the total percent abundance of all organisms for each day adds up to 100%, 
#modify the percent abundance for "Other" for each day in "SRD129_Nasal_FluControlNoDNEGGenusPercentAbundanceAbove1percent.csv" in a spreadsheet editor 
#and save as "SRD129_Nasal_FluControlGenusPercentAbundanceAbove1PercentAddTo100FINAL.csv"

#Create nasal genera plot
nasalgen.2 = read.csv("SRD129_Nasal_FluControlGenusPercentAbundanceAbove1PercentAddTo100FINAL.csv", header = TRUE)
levels(nasalgen.2$Day)  # "D0"  "D1"  "D10" "D14" "D21" "D28" "D3"  "D36" "D42" "D7" 
nasalgen.2$Day = factor(nasalgen.2$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))

#Nasalgen.2 abundance plot for each day, flu and control only, no DNEG12/6, more than 1% genera
(nasalgen.2plot <- ggplot(data=nasalgen.2, aes(x=Treatment, y=Percent.abundance, fill=Genus)) +
    geom_bar(stat = 'identity') +
    #geom_bar(stat= 'identity', colour='black') +
    #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
    facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank()) +
    scale_fill_igv(name = "Genus") +
    theme(legend.direction = "vertical") +
    theme(legend.text = element_text(face = 'italic')))

ggsave("FluControl_NasalGenera.tiff", plot=nasalgen.2plot, width = 15, height = 7, dpi = 500, units =c("in"))
