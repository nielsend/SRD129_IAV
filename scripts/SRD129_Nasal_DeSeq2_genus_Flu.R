#####################################################################################################
#SRD129 Nasal DESeq2 - Genus
#Kathy Mou

#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially abundant between IAV and control groups
#GENUS ONLY, no bad samples

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")

#Annotation
#IC = IAV, control

sessionInfo()
#R version 3.6.3 (2020-02-29)

####### PREPARING OBJECTS FOR DESEQ2 ANALYSIS ########

#Load files
otu2 <- import_mothur(mothur_shared_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.shared')
taxo2 <- import_mothur(mothur_constaxonomy_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
taxo2
meta2 <- read.table(file = './data/SRD129metadata.csv', sep = ',', header = TRUE)

#Organize meta file
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
SRD129 <- prune_samples(sample_sums(SRD129) > 2000, SRD129)  #This removes samples that have fewer than 2000 sequences associated with them.
SRD129 <- prune_taxa(taxa_sums(SRD129) > 10, SRD129)        #Removes OTUs that occur less than 10 times globally
tax_table(SRD129) [1:5, 1:6]  #See what's in tax_table first 5 rows, first 6 columns
SRD129.genus <- tax_glom(SRD129, taxrank = "Genus")
#This method merges species that have the same taxonomy at a certain taxanomic rank. 
#Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 



# PRIMARY COMPARISONS TO MAKE #

# NMDS plot showed that dispersion is different between days, so I subsetted by day

# Important comparisons to make (significant changes in beta diversity between treatments): Compare Days 14, 21, 28 IAV and Control

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
#IAV = 10
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
#IAV = 10
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
#IAV = 10
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
#IAV = 10
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
#IAV = 10
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
#IAV = 10
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
#IAV = 9
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

################################################## Day 28 Nasal ############################################################

sample_data(SRD129.genus)
SRD129.D28 <- subset_samples(SRD129.genus, Day == 'D28')
sample_sums(SRD129.D28)
colnames(otu_table(SRD129.D28))
SRD129.D28 <- prune_taxa(taxa_sums(SRD129.D28) > 1, SRD129.D28)
rowSums(SRD129.D28@otu_table)
SRD129.D28.De <- phyloseq_to_deseq2(SRD129.D28, ~ Set)
SRD129.D28.De <- DESeq(SRD129.D28.De, test = "Wald", fitType = "parametric")

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D28_IAV")
#IAV = 5
sum(meta2$Set == "D28_Control")
#Control = 4

#Extract results from a DESeq analysis, organize table
SRD129.D28.De$Set
res.D28.ic = results(SRD129.D28.De, contrast=c("Set","D28_IAV","D28_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D28.ic = res.D28.ic[which(res.D28.ic$padj < .05), ]
sigtab.D28.ic = cbind(as(sigtab.D28.ic, "data.frame"), as(tax_table(SRD129.D28)[rownames(sigtab.D28.ic), ], "matrix"))
format(sigtab.D28.ic$padj, scientific = TRUE)
sigtab.D28.ic$newp <- format(round(sigtab.D28.ic$padj, digits = 3), scientific = TRUE)
sigtab.D28.ic$Treatment <- ifelse(sigtab.D28.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D28.ic
sum.sigtab.D28.ic <- summary(sigtab.D28.ic)
sum.sigtab.D28.ic

#ggplot
deseq.D28.ic <- ggplot(sigtab.D28.ic, aes(x=reorder(rownames(sigtab.D28.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D28.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 28')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(values = c(Control = "#999999", IAV = "#CC0066"))
deseq.D28.ic

#Add OTU and comparisons columns
sigtab.D28.ic
sigtab.D28.ic$OTU <- rownames(sigtab.D28.ic)
sigtab.D28.ic
sigtab.D28.ic$comp <- 'D28_IAVvsControl'

#Create final significant comparisons table
final.sigtab <- rbind(final.sigtab, sigtab.D28.ic)

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
#IAV = 9
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
#IAV = 9
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

#IAV and Control Log2fold Plot Part A
final_ic <- sigtab.D0.ic
final_ic <- rbind(final_ic, sigtab.D1.ic, sigtab.D3.ic, sigtab.D7.ic, sigtab.D10.ic, 
                  sigtab.D14.ic, sigtab.D21.ic, sigtab.D28.ic, sigtab.D36.ic, sigtab.D42.ic)
final_ic$Family_Genus <- paste(final_ic$Family, final_ic$Genus) #create new column with Family_Genus
final_ic$comp
class(final_ic)
final_ic$comp <- factor(final_ic$comp, levels=c("D0_IAVvsControl",
                                                "D1_IAVvsControl", "D3_IAVvsControl", "D7_IAVvsControl",
                                                "D10_IAVvsControl", "D14_IAVvsControl", "D21_IAVvsControl",
                                                "D28_IAVvsControl", "D36_IAVvsControl", "D42_IAVvsControl"))
levels(final_ic$comp)
ic_plot <- ggplot(final_ic, aes(x=Family_Genus, log2FoldChange, fill = comp)) +
  geom_bar(stat='identity') +
  labs(x="Family Genus", y = "Total log2 Fold Change") +
  theme(axis.text.x=element_text(color = 'black', size = 18),
        axis.text.y=element_text(color = 'black', size=15), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20))+ 
  coord_flip() +
  scale_fill_igv(name = "comp") +
  ggtitle('Differentially Abundant Nasal Families and Genera between Influenza A Virus and Control Groups') + 
  theme(plot.title = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20))
ic_plot <- ic_plot + guides(fill=guide_legend(title="Day and treatment group"))
ic_plot

write.csv(final_ic, file= "SRD129_FinalDiffAbundNasalGenus_IAVcontrol.csv")

#Modified SRD129_FinalDiffAbundNasalGenus_IAVcontrol.csv in a spreadsheet editor by removing all genera except for Actinobacillus, Moraxella, Neisseriaceae_unclassified, Prevotellaceae_NK3B31_group,
#Prevotellaceae_unclassified, Staphylococcus, Streptococcus.
#Saved modified file as SRD129_FinalDiffAbundNasalGenus_IAVcontrolSelectList.csv

#IAV and control Log2fold Plot Part B
ic <- read.csv('SRD129_FinalDiffAbundNasalGenus_IAVcontrolSelectList.csv', header = TRUE, sep = ",")
head(ic[,1:10])
colnames(ic)
ic$DayComp <- sub('_[A-Za-z]+', '\\2', ic$comp)
unique(ic$DayComp)
ic$Day <- sub('_[A-Za-z]+', '\\2', ic$DayComp)
unique(ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D28"
ic$Day = factor(ic$Day, levels=c("D3","D7", "D10","D14", "D21", "D28"))
levels(ic$Day) #"D3"  "D7"  "D10" "D14" "D21" "D28"
(ic_logfoldplot <- ggplot(data=ic, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 3, scales = "free") + ylab('log2-fold change') +
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
ic_logfoldplot <- ic_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
ic_logfoldplot

#Save 'ic_logfoldplot' as a .tiff for publication, 500dpi
ggsave("Figure_5.tiff", plot=ic_logfoldplot, width = 15, height = 7, dpi = 500, units =c("in"))
