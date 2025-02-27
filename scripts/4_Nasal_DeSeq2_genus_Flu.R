########################################
#SRD129 16S - Differential Abundance of Genera in Nasal Microbiota using DESeq2



#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially abundant between treatment groups and control group

#Files needed:
#Mothur shared file: SRD129Flu.outsingletons.abund.opti_mcc.shared
#Mothur constaxonomy file: SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#Metadata: SRD129metadata.csv

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library(cowplot)

#Annotation
#IC = IAV, control

#R version 4.0.2 (2020-06-22)

#######################################################################

#Preparing objects for DESeq2: load files
otu2 <- import_mothur(mothur_shared_file = './SRD129Flu.outsingletons.abund.opti_mcc.shared')
taxo2 <- import_mothur(mothur_constaxonomy_file = './SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
meta2 <- read.table(file = './SRD129metadata.csv', sep = ',', header = TRUE)

#Organize meta file
colnames(meta2) <- c("Sample", "Pig", "Tissue", "Day", "Treatment")
rownames(meta2) <- meta2$Sample
meta2 <- meta2[,-1] #Remove Sample column
meta2$Set <- paste(meta2$Day, meta2$Treatment, sep = '_') #Create 'Set' column that combines Day and Treatment

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta2 <- sample_data(meta2)
SRD129 <- phyloseq(otu2, taxo2)
SRD129 <- merge_phyloseq(SRD129, phy_meta2)   #Combines the metadata with this phyloseq object
colnames(tax_table(SRD129))
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
SRD129

#Prune and subset by genus rank
SRD129 <- prune_samples(sample_sums(SRD129) > 4276, SRD129)  #This removes samples that have fewer than 2000 sequences associated with them.
SRD129 <- prune_taxa(taxa_sums(SRD129) > 10, SRD129)        #Removes OTUs that occur less than 10 times globally
tax_table(SRD129) [1:5, 1:6]  #See what's in tax_table first 5 rows, first 6 columns
SRD129.genus <- tax_glom(SRD129, taxrank = "Genus")
#This method merges species that have the same taxonomy at a certain taxanomic rank.
#Its approach is analogous to tip_glom, but uses categorical data instead of a tree.



# PRIMARY COMPARISONS TO MAKE #

# NMDS plot showed that dispersion is different between days, so I subsetted by day

# Important comparisons to make (significant changes in beta diversity between treatments): Compare Days 14 and 21 IAV and Control

# Other comparisons to make (no significant changes in beta diversity between treatments): Compare Days 0, 1, 3, 7, 10, 36, 42 IAV and Control


################################################## Day 0 Nasal #########################################################
  
sample_data(SRD129.genus)
SRD129.D0 <- subset_samples(SRD129.genus, Day == 'D0')
sample_sums(SRD129.D0)
colnames(otu_table(SRD129.D0)) #Check on all the sample names
SRD129.D0 <- prune_taxa(taxa_sums(SRD129.D0) > 1, SRD129.D0)
#If taxa_sums is >1, then it will print that out in SRD129.D0 object and not include anything with <1.
rowSums(SRD129.D0@otu_table)
SRD129.D0.De <- phyloseq_to_deseq2(SRD129.D0, ~ Set)
# ~Set: whatever you want to group data by and use to designate ellipses with
SRD129.D0.De <- DESeq(SRD129.D0.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using "meta2" dataframe):
sum(meta2$Set == "D0_IAV")
#IAV = 7
sum(meta2$Set == "D0_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D0.De$Set
res.D0.ic = results(SRD129.D0.De, contrast=c("Set","D0_IAV","D0_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0.ic = res.D0.ic[which(res.D0.ic$padj < .05), ]
sigtab.D0.ic = cbind(as(sigtab.D0.ic, "data.frame"), as(tax_table(SRD129.D0)[rownames(sigtab.D0.ic), ], "matrix"))
format(sigtab.D0.ic$padj, scientific = TRUE)
sigtab.D0.ic$newp <- format(round(sigtab.D0.ic$padj, digits = 3), scientific = TRUE)
sigtab.D0.ic$Treatment <- ifelse(sigtab.D0.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D0.ic
sum.sigtab.D0.ic <- summary(sigtab.D0.ic)
sum.sigtab.D0.ic

#ggplot
deseq.D0.ic <- ggplot(sigtab.D0.ic, aes(x=reorder(rownames(sigtab.D0.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control="#999999", IAV = "#CC0066"))
deseq.D0.ic

#Add OTU and comparisons columns
sigtab.D0.ic
sigtab.D0.ic$OTU <- rownames(sigtab.D0.ic)
sigtab.D0.ic
sigtab.D0.ic$comp <- 'D0_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(sigtab.D0.ic)

################################################## Day 1 Nasal ######################################################################

sample_data(SRD129.genus)
SRD129.D1 <- subset_samples(SRD129.genus, Day == 'D1')
sample_sums(SRD129.D1)
colnames(otu_table(SRD129.D1)) #check on all the sample names
SRD129.D1 <- prune_taxa(taxa_sums(SRD129.D1) > 1, SRD129.D1)
rowSums(SRD129.D1@otu_table)
SRD129.D1.De <- phyloseq_to_deseq2(SRD129.D1, ~ Set)
SRD129.D1.De <- DESeq(SRD129.D1.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D1_IAV")
#IAV = 7
sum(meta2$Set == "D1_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D1.De$Set
res.D1.ic = results(SRD129.D1.De, contrast=c("Set","D1_IAV","D1_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D1.ic = res.D1.ic[which(res.D1.ic$padj < .05), ]
sigtab.D1.ic = cbind(as(sigtab.D1.ic, "data.frame"), as(tax_table(SRD129.D1)[rownames(sigtab.D1.ic), ], "matrix"))
format(sigtab.D1.ic$padj, scientific = TRUE)
sigtab.D1.ic$newp <- format(round(sigtab.D1.ic$padj, digits = 3), scientific = TRUE)
sigtab.D1.ic$Treatment <- ifelse(sigtab.D1.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D1.ic
sum.sigtab.D1.ic <- summary(sigtab.D1.ic)
sum.sigtab.D1.ic

#ggplot
deseq.D1.ic <- ggplot(sigtab.D1.ic, aes(x=reorder(rownames(sigtab.D1.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D1.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 1')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D1.ic

#Add OTU and comparisons columns
sigtab.D1.ic
sigtab.D1.ic$OTU <- rownames(sigtab.D1.ic)
sigtab.D1.ic
sigtab.D1.ic$comp <- 'D1_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D1.ic)

################################################## Day 3 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D3 <- subset_samples(SRD129.genus, Day == 'D3')
sample_sums(SRD129.D3)
colnames(otu_table(SRD129.D3)) #check on all the sample names
SRD129.D3 <- prune_taxa(taxa_sums(SRD129.D3) > 1, SRD129.D3)
rowSums(SRD129.D3@otu_table)
SRD129.D3.De <- phyloseq_to_deseq2(SRD129.D3, ~ Set)
SRD129.D3.De <- DESeq(SRD129.D3.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D3_IAV")
#IAV = 7
sum(meta2$Set == "D3_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D3.De$Set
res.D3.ic = results(SRD129.D3.De, contrast=c("Set","D3_IAV","D3_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D3.ic = res.D3.ic[which(res.D3.ic$padj < .05), ]
sigtab.D3.ic = cbind(as(sigtab.D3.ic, "data.frame"), as(tax_table(SRD129.D3)[rownames(sigtab.D3.ic), ], "matrix"))
format(sigtab.D3.ic$padj, scientific = TRUE)
sigtab.D3.ic$newp <- format(round(sigtab.D3.ic$padj, digits = 3), scientific = TRUE)
sigtab.D3.ic$Treatment <- ifelse(sigtab.D3.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D3.ic
sum.sigtab.D3.ic <- summary(sigtab.D3.ic)
sum.sigtab.D3.ic

#ggplot
deseq.D3.ic <- ggplot(sigtab.D3.ic, aes(x=reorder(rownames(sigtab.D3.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D3.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 3')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D3.ic

#Add OTU and comparisons columns
sigtab.D3.ic
sigtab.D3.ic$OTU <- rownames(sigtab.D3.ic)
sigtab.D3.ic
sigtab.D3.ic$comp <- 'D3_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D3.ic)

################################################## Day 7 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D7 <- subset_samples(SRD129.genus, Day == 'D7')
sample_sums(SRD129.D7)
colnames(otu_table(SRD129.D7)) #check on all the sample names
SRD129.D7 <- prune_taxa(taxa_sums(SRD129.D7) > 1, SRD129.D7)
rowSums(SRD129.D7@otu_table)
SRD129.D7.De <- phyloseq_to_deseq2(SRD129.D7, ~ Set)
SRD129.D7.De <- DESeq(SRD129.D7.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D7_IAV")
#IAV = 7
sum(meta2$Set == "D7_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D7.De$Set
res.D7.ic = results(SRD129.D7.De, contrast=c("Set","D7_IAV","D7_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D7.ic = res.D7.ic[which(res.D7.ic$padj < .05), ]
sigtab.D7.ic = cbind(as(sigtab.D7.ic, "data.frame"), as(tax_table(SRD129.D7)[rownames(sigtab.D7.ic), ], "matrix"))
format(sigtab.D7.ic$padj, scientific = TRUE)
sigtab.D7.ic$newp <- format(round(sigtab.D7.ic$padj, digits = 3), scientific = TRUE)
sigtab.D7.ic$Treatment <- ifelse(sigtab.D7.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D7.ic
sum.sigtab.D7.ic <- summary(sigtab.D7.ic)
sum.sigtab.D7.ic

#ggplot
deseq.D7.ic <- ggplot(sigtab.D7.ic, aes(x=reorder(rownames(sigtab.D7.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D7.ic

#Add OTU and comparisons columns
sigtab.D7.ic
sigtab.D7.ic$OTU <- rownames(sigtab.D7.ic)
sigtab.D7.ic
sigtab.D7.ic$comp <- 'D7_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D7.ic)

################################################## Day 10 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D10 <- subset_samples(SRD129.genus, Day == 'D10')
sample_sums(SRD129.D10)
colnames(otu_table(SRD129.D10)) #check on all the sample names
SRD129.D10 <- prune_taxa(taxa_sums(SRD129.D10) > 1, SRD129.D10)
rowSums(SRD129.D10@otu_table)
SRD129.D10.De <- phyloseq_to_deseq2(SRD129.D10, ~ Set)
SRD129.D10.De <- DESeq(SRD129.D10.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D10_IAV")
#IAV = 7
sum(meta2$Set == "D10_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D10.De$Set
res.D10.ic = results(SRD129.D10.De, contrast=c("Set","D10_IAV","D10_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D10.ic = res.D10.ic[which(res.D10.ic$padj < .05), ]
sigtab.D10.ic = cbind(as(sigtab.D10.ic, "data.frame"), as(tax_table(SRD129.D10)[rownames(sigtab.D10.ic), ], "matrix"))
format(sigtab.D10.ic$padj, scientific = TRUE)
sigtab.D10.ic$newp <- format(round(sigtab.D10.ic$padj, digits = 3), scientific = TRUE)
sigtab.D10.ic$Treatment <- ifelse(sigtab.D10.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D10.ic
sum.sigtab.D10.ic <- summary(sigtab.D10.ic)
sum.sigtab.D10.ic

#ggplot
deseq.D10.ic <- ggplot(sigtab.D10.ic, aes(x=reorder(rownames(sigtab.D10.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D10.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 10')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D10.ic

#Add OTU and comparisons columns
sigtab.D10.ic
sigtab.D10.ic$OTU <- rownames(sigtab.D10.ic)
sigtab.D10.ic
sigtab.D10.ic$comp <- 'D10_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D10.ic)

################################################## Day 14 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D14 <- subset_samples(SRD129.genus, Day == 'D14')
sample_sums(SRD129.D14)
colnames(otu_table(SRD129.D14)) #check on all the sample names
SRD129.D14 <- prune_taxa(taxa_sums(SRD129.D14) > 1, SRD129.D14)
rowSums(SRD129.D14@otu_table)
SRD129.D14.De <- phyloseq_to_deseq2(SRD129.D14, ~ Set)
SRD129.D14.De <- DESeq(SRD129.D14.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D14_IAV")
#IAV = 7
sum(meta2$Set == "D14_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D14.De$Set
res.D14.ic = results(SRD129.D14.De, contrast=c("Set","D14_IAV","D14_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D14.ic = res.D14.ic[which(res.D14.ic$padj < .05), ]
sigtab.D14.ic = cbind(as(sigtab.D14.ic, "data.frame"), as(tax_table(SRD129.D14)[rownames(sigtab.D14.ic), ], "matrix"))
format(sigtab.D14.ic$padj, scientific = TRUE)
sigtab.D14.ic$newp <- format(round(sigtab.D14.ic$padj, digits = 3), scientific = TRUE)
sigtab.D14.ic$Treatment <- ifelse(sigtab.D14.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D14.ic
sum.sigtab.D14.ic <- summary(sigtab.D14.ic)
sum.sigtab.D14.ic

#ggplot
deseq.D14.ic <- ggplot(sigtab.D14.ic, aes(x=reorder(rownames(sigtab.D14.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 14')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D14.ic

#Add OTU and comparisons columns
sigtab.D14.ic
sigtab.D14.ic$OTU <- rownames(sigtab.D14.ic)
sigtab.D14.ic
sigtab.D14.ic$comp <- 'D14_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D14.ic)

################################################## Day 21 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D21 <- subset_samples(SRD129.genus, Day == 'D21')
sample_sums(SRD129.D21)
colnames(otu_table(SRD129.D21))
SRD129.D21 <- prune_taxa(taxa_sums(SRD129.D21) > 1, SRD129.D21)
rowSums(SRD129.D21@otu_table)
SRD129.D21.De <- phyloseq_to_deseq2(SRD129.D21, ~ Set)
SRD129.D21.De <- DESeq(SRD129.D21.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D21_IAV")
#IAV = 7
sum(meta2$Set == "D21_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D21.De$Set
res.D21.ic = results(SRD129.D21.De, contrast=c("Set","D21_IAV","D21_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.ic = res.D21.ic[which(res.D21.ic$padj < .05), ]
sigtab.D21.ic = cbind(as(sigtab.D21.ic, "data.frame"), as(tax_table(SRD129.D21)[rownames(sigtab.D21.ic), ], "matrix"))
format(sigtab.D21.ic$padj, scientific = TRUE)
sigtab.D21.ic$newp <- format(round(sigtab.D21.ic$padj, digits = 3), scientific = TRUE)
sigtab.D21.ic$Treatment <- ifelse(sigtab.D21.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D21.ic
sum.sigtab.D21.ic <- summary(sigtab.D21.ic)
sum.sigtab.D21.ic

#ggplot
deseq.D21.ic <- ggplot(sigtab.D21.ic, aes(x=reorder(rownames(sigtab.D21.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 21')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D21.ic

#Add OTU and comparisons columns
sigtab.D21.ic
sigtab.D21.ic$OTU <- rownames(sigtab.D21.ic)
sigtab.D21.ic
sigtab.D21.ic$comp <- 'D21_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D21.ic)


################################################## Day 36 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D36 <- subset_samples(SRD129.genus, Day == 'D36')
sample_sums(SRD129.D36)
colnames(otu_table(SRD129.D36))
SRD129.D36 <- prune_taxa(taxa_sums(SRD129.D36) > 1, SRD129.D36)
rowSums(SRD129.D36@otu_table)
SRD129.D36.De <- phyloseq_to_deseq2(SRD129.D36, ~ Set)
SRD129.D36.De <- DESeq(SRD129.D36.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D36_IAV")
#IAV = 7
sum(meta2$Set == "D36_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D36.De$Set
res.D36.ic = results(SRD129.D36.De, contrast=c("Set","D36_IAV","D36_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D36.ic = res.D36.ic[which(res.D36.ic$padj < .05), ]
sigtab.D36.ic = cbind(as(sigtab.D36.ic, "data.frame"), as(tax_table(SRD129.D36)[rownames(sigtab.D36.ic), ], "matrix"))
format(sigtab.D36.ic$padj, scientific = TRUE)
sigtab.D36.ic$newp <- format(round(sigtab.D36.ic$padj, digits = 3), scientific = TRUE)
sigtab.D36.ic$Treatment <- ifelse(sigtab.D36.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D36.ic
sum.sigtab.D36.ic <- summary(sigtab.D36.ic)
sum.sigtab.D36.ic

#ggplot
deseq.D36.ic <- ggplot(sigtab.D36.ic, aes(x=reorder(rownames(sigtab.D36.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D36.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 36')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D36.ic

#Add OTU and comparisons columns
sigtab.D36.ic
sigtab.D36.ic$OTU <- rownames(sigtab.D36.ic)
sigtab.D36.ic
sigtab.D36.ic$comp <- 'D36_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D36.ic)

################################################## Day 42 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D42 <- subset_samples(SRD129.genus, Day == 'D42')
sample_sums(SRD129.D42)
colnames(otu_table(SRD129.D42))
SRD129.D42 <- prune_taxa(taxa_sums(SRD129.D42) > 1, SRD129.D42)
rowSums(SRD129.D42@otu_table)
SRD129.D42.De <- phyloseq_to_deseq2(SRD129.D42, ~ Set)
SRD129.D42.De <- DESeq(SRD129.D42.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe):
sum(meta2$Set == "D42_IAV")
#IAV = 7
sum(meta2$Set == "D42_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D42.De$Set
res.D42.ic = results(SRD129.D42.De, contrast=c("Set","D42_IAV","D42_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D42.ic = res.D42.ic[which(res.D42.ic$padj < .05), ]
sigtab.D42.ic = cbind(as(sigtab.D42.ic, "data.frame"), as(tax_table(SRD129.D42)[rownames(sigtab.D42.ic), ], "matrix"))
format(sigtab.D42.ic$padj, scientific = TRUE)
sigtab.D42.ic$newp <- format(round(sigtab.D42.ic$padj, digits = 3), scientific = TRUE)
sigtab.D42.ic$Treatment <- ifelse(sigtab.D42.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D42.ic
sum.sigtab.D42.ic <- summary(sigtab.D42.ic)
sum.sigtab.D42.ic

#ggplot
deseq.D42.ic <- ggplot(sigtab.D42.ic, aes(x=reorder(rownames(sigtab.D42.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D42.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 42')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D42.ic

#Add OTU and comparisons columns
sigtab.D42.ic
sigtab.D42.ic$OTU <- rownames(sigtab.D42.ic)
sigtab.D42.ic
sigtab.D42.ic$comp <- 'D42_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D42.ic)

################################################## write csv ########################################
write.csv(final.sigtab, file= "SRD129_FinalDiffAbundGenus.csv")


######### Plots of Differentially Abundant Nasal Genera Combined for Each Pairwise Comparison of IAV and Control Groups########
library("ggsci")

#IAV and Control Log2fold Plot
#create new function to use in filter below
`%notin%` <- function(x,y) !(x %in% y) 

#Modified SRD129_FinalDiffAbundGenus.csv by removing all genera that showed difference on day 0
final.sigtab2 <- final.sigtab
sigtab.D0.ic1 <- as.matrix(sigtab.D0.ic)
sigtab.D0.ic1OTU <- row.names(sigtab.D0.ic1)
ic <- final.sigtab2 %>% filter(OTU %notin% sigtab.D0.ic1OTU) #57 OTU dif on day 0

#Finder No of unique genera not dif. in day 0
uniqueGenera <- unique(ic$Genus)
uniqueOTU <- unique(ic$OTU)
View(final.sigtab2)



write.csv(ic, './SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDay0Genera.csv')

#Saved modified file as SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDay0Genera.csv
ic <- read.csv('./SRD129_FinalDiffAbundNasalGenus_IAVcontrol_NoDay0Genera.csv', header = TRUE, sep = ",")
View(ic)
head(ic[,1:10])
colnames(ic)
ic$DayComp <- sub('_[A-Za-z]+', '\\2', ic$comp)
unique(ic$DayComp)
ic$Day <- sub('_[A-Za-z]+', '\\2', ic$DayComp)
unique(ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
ic$Day = factor(ic$Day, levels=c("D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))
levels(ic$Day) #"D1" "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
(ic_logfoldplot <- ggplot(data=ic, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15)) +
    scale_fill_manual(values = c(Control = "#F8766D", IAV = "#00BFC4")))
ic_logfoldplot <- ic_logfoldplot + theme(strip.text = element_text(size= 11, face='italic'))
ic_logfoldplot

#Save 'ic_logfoldplot' as a .tiff for publication, 500dpi
ggsave("Supplemental_Figure_2.tiff", plot=ic_logfoldplot, width = 20, height = 40, dpi = 500, units =c("in"))

library("dplyr")
#Filter out unknown Genera or non-clinically relevant Genera of a single day
ic2 <- ic
ic2$Genus <- as.character(ic2$Genus)
ic3 <- ic4 %>% subset(`Genus` != "vadinBE97_ge" & `Genus` !=  "uncultured_ge" & Genus != "Mitsuokella"& 
                        Genus != "Butyrivibrio" & Genus != "Corynebacteriaceae_unclassified" & 
                        Genus != "p-1088-a5_gut_group" & Genus != "WCHB1-41_ge" & Genus != "Selenomonas_3" &
                        Genus != "Ruminococcaceae_UCG-004" & Genus != "Puniceicoccaceae_unclassified" & 
                        Genus != "Micrococcales_unclassified" & Genus != "Arcanobacterium"& Genus != "Angelakisella"& 
                        Genus != "alphaproteobacteria_unclassified" & Genus != "Acetitomaculum" &  Genus != "Aerosphaera" &
                        Genus != "Anaerostipes" & Genus != "horsej-a03" & Genus != "uncultured" & 
                        Genus != "Firmicutes_unclassified" & Genus != "Family_XIII_AD3011_group" & Genus != "Campylobacter" &
                        Genus != "Bacteria_unclassified" & Genus != "Kurthia" & Genus != "Lachnoclostridium_12" &
                        Genus != "Lachnospiraceae_FCS020_group" & Genus != "Lachnospiraceae_NK4B4_group" &
                        Genus != "Prevotellaceae_UCG-004" & Genus != "Alphaproteobacteria_unclassified" & 
                        Genus != "Clostridiaceae_1-unclassified" & Genus != "Bradymonadales_ge" & 
                        Genus != "CAG-56" & Genus != "Cellulosilyticum" & Genus != "Clostridiaceae_1_unclassified" &
                        Genus != "Actinobacillus" & Genus != "Dyadobacter" & Genus != "Enhydrobacter" & Genus != "Erysipelotrichaceae_UCG-009" & 
                        Genus != "Fastidiosipila" & Genus != "Filifactor" & Genus != "Gemella" & Genus != "Izimaplasmatales_ge" &
                        Genus != "Lachnospiraceae_UCG-004" & Genus != "Mogibacterium" & Genus != "Mollicutes_RF39_ge" & 
                        Genus != "Moraxellaceae_unclassified" & Genus != "Myroides" & Genus != "Negativibacillus" & 
                        Genus != "Neisseriaceae_unclassified" & Genus != "Oribacterium" & Genus != "Parabacteroides" & Genus != "Parvimonas" &
                        Genus != "Pedobacter" & Genus != "Peptococcaceae_unclassified" & Genus != "Peptococcus" & Genus != "Peptoniphilus" &
                        Genus != "PPorphyromonas" & Genus != "Prevotellaceae_unclassified" & Genus != "Rhizobiaceae_unclassified" & Genus != "Ruminiclostridium_6" &
                        Genus != "Ruminiclostridium_9" & Genus != "Ruminococcaceae_ge" & Genus != "Ruminococcus_2" & Genus != "Veillonella" &
                        Genus != "Bacteroides" & Genus != "Bacteroida_unclassified" & Genus != "Allofustis" & Genus != "Aerococcaceae_unclassified"&
                        Genus != "Bergeyella" &  Genus != "Brachybacterium" & Genus != "Clostridales_vadinBB60_group_ge"&
                        Genus != "Deinococcus" & Genus != "Dorea" & Genus != "Fibrobacter" & Genus != "Lachnospiraceae_ge" & 
                        Genus != "Leptotrichia" & Genus !="Porphyomonas" & Genus != "Prevotella_1" & Genus != "Ruminococcaceae_UCG-010" &
                        Genus != "Spaerochaeta" & Genus != "Ruminococcus_1" & Genus != "Christensenellaceae_R-7_group" & Genus != "Clostridiales_vadinBB60_group_ge" &
                        Genus != "Holdemanella" & Genus != "Porphyromonas" & Genus != "Ruminococcaceae_UCG-008" & Genus != "Ruminococcaceae_unclassified" & Genus != "Sphaerochaeta" & 
                        Genus != "Bacteroidia_unclassified" & Genus != "Candidatus_Soleaferrea" & Genus != "" & Genus != "" &
                        Genus != "Weissella" & Genus != "Vagococcus" & Genus != "Subdoligranulum" & Genus != "Streptococcus" & 
                        Genus != "Sphingobacterium" & Genus != "Pseudomonas" & Genus != "Prevotellaceae_UCG-003" & 
                        Genus != "Prevotella_2" & Genus != "Massilla" & Genus != "Fusobacteriaceae_unclassified" & Genus != "Faecalibacterium" & 
                        Genus != "Coprococcus_1" & Genus != "Collinsetia" & Genus != "Cateribacterium" & Genus != "Brevundimonas" & 
                        Genus != "Blautia" & Genus != "Acidaminococcus" & Genus != "Actinomyces" & Genus != "Catenisphaera" &
                        Genus != "Collinsella" & Genus != "Enterohabidus" & Genus != "Lechnospiraceae_NKA320_group" & 
                        Genus != "Massilia" & Genus != "Pseudoclavibacter" & Genus != "Ruminococcacae_UCG-014" & Genus != "Succiniclasticum" & 
                        Genus != "Wohlfahrtiimonas" & Genus != "Actinomyces" & Genus != "Burkholderiaceae_unclassified" & Genus != "Enterorhabdus" & 
                        Genus != "Lachnospiraceae_NK3A20_group" & Genus != "Ruminococcaceae_UCG-014" & Genus != "" & Genus != "")


ic4 <- ic3 %>% select(Day, log2FoldChange, Treatment, Genus)
(ic3_logfoldplot <- ggplot(data=ic3, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15)) +
    scale_fill_manual(values = c(Control = "#F8766D", IAV = "#00BFC4")))
ic3_logfoldplot <- ic3_logfoldplot + theme(strip.text = element_text(size= 12, face='italic'))
ic3_logfoldplot
gc()
ggsave("Supplemental_Figure_2-4.tiff", plot=ic3_logfoldplot, width = 25, height = 30, dpi = 500, units =c("in"))
unique(ic$Genus)


# Examine the things that WERE different on Day 0
ic.day0 <- final.sigtab2 %>% filter(OTU %in% sigtab.D0.ic1OTU) #55 OTU
unique(ic.day0$Genus)

ic.day0$DayComp <- sub('_[A-Za-z]+', '\\2', ic.day0$comp)
unique(ic.day0$DayComp)
ic.day0$Day <- sub('_[A-Za-z]+', '\\2', ic.day0$DayComp)
unique(ic.day0$Day) #"D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
ic.day0$Day = factor(ic.day0$Day, levels=c("D0", "D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))
levels(ic.day0$Day) #"D1" "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
(ic.day0_logfoldplot <- ggplot(data=ic.day0, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15)) +
    scale_fill_manual(values = c(Control = "#F8766D", IAV = "#00BFC4")))
gc()
ic.day0_logfoldplot <- ic.day0_logfoldplot + theme(strip.text = element_text(size= 11, face='italic'))
ic.day0_logfoldplot

#Save 'ic_logfoldplot' as a .tiff for publication, 500dpi
ggsave("Supplemental_Figure_2_day0dif.tiff", plot=ic.day0_logfoldplot, width = 20, height = 25, dpi = 500, units =c("in"))

library(tidyverse)



#Examine Respiratory pathogens of interest
RPI.ic <- final.sigtab2
sample_sums(RPI.ic)
colnames(otu_table(RPI.ic)) #check on all the sample names
RPI.ic <- prune_taxa(taxa_sums(RPI.ic) > 1, RPI.ic)
rowSums(RPI.ic@otu_table)
RPI.ic <- phyloseq_to_deseq2(RPI.ic, ~ Set)
RPI.ic <- DESeq(RPI.ic, test = "Wald", fitType = "parametric")


RPI.ic$DayComp <- sub('_[A-Za-z]+', '\\2', RPI.ic$comp)
unique(RPI.ic$DayComp)
RPI.ic$Day <- sub('_[A-Za-z]+', '\\2', RPI.ic$DayComp)
unique(RPI.ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
RPI.ic$Day = factor(RPI.ic$Day, levels=c("D0", "D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))
levels(RPI.ic$Day) #"D1" "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"

RPI <- RPI.ic %>% subset(`Genus` == "Actinobacillus" | `Genus` == "Mycoplasma" | `Genus` == "Pasteurellaceae_unclassified" | `Genus` == "Streptococcus" |  `Genus` == "Trueperella" | `Genus` == "Salmonella")
levels(RPI$Day)
write.csv(RPI, "./RPI.csv") #add Day 1 & 14 for graph
RPI <- read.csv("./RPI_dayAdd.csv")
colnames(RPI) <- c("log2FoldChange", "Phylum", "Class", "Order", "Family", "Genus", "DayComp", "Day", "Treatment")
RPI$Day  = factor(RPI$Day, levels=c("D0", "D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))

(RPI_logfoldplot <- ggplot(data=RPI, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5) + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15)) +
    scale_fill_manual(values = c(Control = "#F8766D", IAV = "#00BFC4")))
RPI_logfoldplot <- RPI_logfoldplot + theme(strip.text = element_text(size= 12, face='italic'))
RPI_logfoldplot
gc()
ggsave("RPI_fig.tiff", plot=RPI_logfoldplot, width = 20, height = 10, dpi = 500, units =c("in"))

# Merge plots
log2FoldPlots <- plot_grid(ic3_logfoldplot + theme(legend.position = "none"), RPI_logfoldplot + theme(legend.position = "bottom"), labels = "AUTO", rel_widths = c(50,50), rel_heights= c(80, 20), scale=0.9, align=c("v"), axis=c("t"), ncol=1)
log2FoldPlots
?plot_grid
ggsave("log2FoldPlots.tiff", plot=log2FoldPlots, width = 20, height = 20, dpi = 500, units =c("in"))



