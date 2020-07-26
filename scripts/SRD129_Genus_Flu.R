############################################################
#SRD129 IAV, Control total genera
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs
#Only analyzed data without bad samples

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

source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)


#Import files
otu <- import_mothur(mothur_shared_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('./data/SRD129Flu.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('./data/SRD129metadata.csv', header = TRUE, sep = ",")
head(meta)



###########
#How to make fobar2.gather file?
#Check out phyloseq R scripts

colnames(meta)[1] <- 'group' 
#Rename first column of "meta" as "group" temporarily. Will use "group" to set as rownames later and remove the "group" column
meta$Day<- gsub("D", "", meta$Day) #Remove "D"
meta$group <- as.character(meta$group)
head(meta)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
head(phy_meta)
phy_meta <- phy_meta[,-1]
head(phy_meta)

#Create phyloseq-class objects with "otu" and "taxo"
FS1b <- phyloseq(otu, taxo)
FS1b <- merge_phyloseq(FS1b, phy_meta)  #This combines the 'phy_meta' metadata with 'FS1b' phyloseq object
colnames(tax_table(FS1b)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_sums(FS1b) #Calculate the sum of all OTUs for each sample. Almost all samples have 2000 sequences
FS1b <- prune_taxa(taxa_sums(FS1b) > 2, FS1b)  #Removes OTUs that occur less than 2 times globally
FS1b.genus <- tax_glom(FS1b, 'Genus')
phyla_tab <- as.data.frame(t(FS1b.genus@otu_table)) #Transpose 'Fs1b.genus' by "otu_table"
head(phyla_tab)
FS1b.genus@tax_table[,6]
colnames(phyla_tab) <- FS1b.genus@tax_table[,6] #Replace column names in phyla_tab from Otuxxxx with Genus names
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #Calculate the proportion of specific phyla per phyla column in 'phyla_tab'
head(phyla_tab2)
phyla_tab2$group <- rownames(phyla_tab2) #Create new column called "group" in 'phyla_tab2' containing rownames
head(phyla_tab2)
fobar <- merge(meta, phyla_tab2, by = 'group') #Merge 'meta' with 'phyla_tab2' by "group"
head(fobar)
fobar.gather <- fobar %>% gather(Genus, value, -(group:Treatment))  #This converts 'fobar' to long-form dataframe. This is handy for using ggplot faceting functions, check out tidyverse tutorials
#This also created new columns "Genus", "value"; it added columns "group" through "Treatment" before "Genus" and "value"
head(fobar.gather)

#Check to see if there is an extra "group" column. If so, run the next set of commands (up to "head(fobar2)") and remove appropriate column
which(colnames(phyla_tab2)=="group") #Results say column 237 is "group" column
phyla_tab3 <- phyla_tab2[,-237] #Drop the 237th column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #Keep the columns that have greater than 0.1 value
phyla_tab4$group <- rownames(phyla_tab4) #Rename rownames as "group"
fobar2 <- merge(meta, phyla_tab4, by = 'group')
head(fobar2)

#To see how many TT are in meta$Tissue: 
length(which(meta$Tissue== "TT")) 
#Output:
#35

fobar2.gather$day <- NULL

#Reorder days 0-14 in 'fobar2.gather' plot
levels(sample_data(fobar2.gather)$Day)
fobar2.gather$Day <- factor(fobar2.gather$Day, levels=c("D0", "D4", "D7", "D11", "D14"))
head(fobar2.gather$Day)

#Create "All" column with "Day", "Treatment" and "Tissue" in 'fobar2.gather'
fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, fobar2.gather$Tissue, sep = '_')

#Count the number of unique items in 'fobar2.gather'. We're interested in the total unique number of genera
fobar2.gather %>% summarise_each(funs(n_distinct)) #90 total unique genera
fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/90))*100) #90 refers to number of Genera

#Subset each "All" group from fobar2.gather and save as individual csv files.
#Two examples:
D0_Control_NW.genus <- subset(fobar2.gather, All=="D0_Control_NW")
write.csv(D0_Control_NW.genus, file="D0_Control_NW.genus.csv")

D4_Infeed_TT.genus <- subset(fobar2.gather, All=="D4_Infeed_TT.genus")
write.csv(D4_Infeed_TT.genus, file="D4_Infeed_TT.genus")

#############

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
