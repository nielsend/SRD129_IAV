#########################################
#SRD129 16S phyloseq prep, adonis, PCoA each tissue- flu and control
#Kathy Mou

#NOTES: 
#This code creates phyloseq objects ("All" or no "All" column) for beta diversity statistics 
#of nasal and tonsil tissue samples separately
## a couple adonis functions and results that you should include in paper (treatment,
#day, day&treatment effects on variation)
#PCoA plots (day, treatment, tissue if applicable, treatment x day) to see general trends for each tissue

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory on desktop
#Mac
setwd("~/Desktop/SRD129/SRD129_2000singletons")

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)

#Load image file
load("SRD129_phyloseq_each_tissue_nobad2.RData")
load("SRD129_phyloseq_Flu.RData")

#Save image file
save.image(file="SRD129_phyloseq_each_tissue_nobad2.RData")
save.image(file="SRD129_phyloseq_Flu.RData")
#########################################

#Read files
otu <- read.csv("SRD129abundsingleton2000OTUtable.csv", row.names=1) #set column 1 as row names
meta <-read.csv("SRD129metadata.csv")
#meta <-read.csv("SRD129metadata06292018.csv")

#Remove the following samples using a text editor from SRD129abundsingleton2000OTUtable.csv 
#and SRD129metadata.csv because these samples were accidentally contaminated during 
#DNA extraction process: 
#SRD129_P150_N_D28 
#SRD129_P151_N_D28 
#SRD129_P152_N_D28 
#SRD129_P153_N_D28 
#SRD129_P154_N_D28 
#SRD129_P156_N_D28 
#SRD129_P163_N_D28 
#SRD129_P164_N_D28 
#SRD129_P165_N_D28 
#SRD129_P166_N_D28 
#SRD129_P167_N_D28 
#SRD129_P168_N_D28 

dim(otu) #1299 187
head(otu[,1:5])
head(otu[,185:187])

dim(meta) #186 5
head(meta[,1:5])

#Remove taxonomy
tax <- otu[,(186:187)] #removed column 186 taxonomy and 187 to include the row names
head(tax)
colnames(tax)[1] <- "delete" #renamed column 1 (formerly 908) as "delete" which will later be deleted
head(tax)

#OTU count data only
otu <- otu[,-187] #remove column 187 taxonomy to have only otu data
head(otu[,180:186]) 
dim(otu) #1299 rows 186 columns

#Transpose to match metadata format
otu.trans <- t(otu) #now rownames are sample names, columns are OTUs
head(otu.trans[,1:5])
head(meta) #row names are numbered, but want sample names as row names
rownames(meta) <- meta$Sample #Set names under "Sample" as rownames
head(meta)
meta <- meta[,-1]
head(meta)

class(meta) #dataframe
class(otu) #dataframe

#Merge otu and meta data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) 
#x=0 means match via rownames from meta; y=0 means match via rownames from otu.trans
head(otu.meta[,1:10])
class(otu.meta) #dataframe
dim(otu.meta) #186 1304

otu.meta2<- cbind(otu.meta) #Make second copy of 'otu.meta' to use to include an "All" column
otu.meta2$All <- with(otu.meta2, paste0(Day, sep=" ", Treatment)) #Combine "Day" and "Treatment" columns into an "All" column
head(otu.meta2) #Check first part of otu.meta2
dim(otu.meta2) #Check dimensions of 'otu.meta2' dataframe
head(otu.meta2[,1300:1305]) #Check the first part of the end of 'otu.meta2' dataframe
head(otu.meta2[,1:10])
rownames(otu.meta2) <- otu.meta2$Row.names
head(otu.meta2[,1:10])
otu.meta2 <- otu.meta2[,-1]
head(otu.meta2[,1:10])
dim(otu.meta2) #186 1304
head(otu.meta2[,1300:1304])
otu.meta2<- otu.meta2[,c(1:4,1304,5:1303)] #Reorder columns to have "All" column after "Treatment" column
head(otu.meta2[,1:10]) #Check the first part of the beginning of 'otu.meta2' dataframe
head(otu.meta2[,1300:1303])
dim(otu.meta2)
write.csv(otu.meta2, file="SRD129abundsingleton2000.otu.meta.csv")


##############################-------NASAL PHYLOSEQ with "All" column--------##########################

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
dim(otu.meta2) #185 1304
otu.meta3 <- otu.meta2[,1:5] #Take columns 1-5 to make 'otu.meta3'
head(otu.meta3)
dim(otu.meta3) #186 5

#Create SAM metada table phyloseq object
SAM = sample_data(otu.meta3, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #186 5
unique(SAM$All)
nrow(SAM[SAM$Day == "D28",]) #9

#OTU only
head(otu.meta2[,1:10])
dim(otu.meta2) #186 1304
head(otu.meta2[,1300:1304])
otu.meta4 <- otu.meta2[,c(6:1304)] #Select "OTU" columns to create 'otu.meta4' dataframe
head(otu.meta4[,1:10])
dim(otu.meta4) #186 1299
otu.meta4.trans <- t(otu.meta4) #Transpose 'otu.meta4' to have OTUs as rownames, sample names as column names
head(otu.meta4.trans[,1:10])
dim(otu.meta4.trans) #1299  186
head(otu.meta4.trans[,180:186])

#Merge 'tax' back into 'otu.meta4.trans' for correct format and taxons
head(tax)
otu.tax <- merge(otu.meta4.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu.tax) #1299  189
head(otu.tax[,185:189])
head(otu.tax[,1:5])
row.names(otu.tax) <- otu.tax[,1] #Set first row as rownames
head(otu.tax[,1:5])
otu.tax <- otu.tax[,-1] #Remove first row, extraneous OTU column
head(otu.tax[,1:5])
dim(otu.tax) #1299 188
head(otu.tax[185:188])

#Split again
dim(otu.tax) #1299 188
head(otu.tax[185:188])
head(otu.tax[1:10])
otu.notax <- otu.tax[,1:186] #Take rows 1-186 to make new dataframe 'otu.notax' (187 is delete column, 188 is taxonomy column)
head(otu.notax[,1:5])
head(otu.notax[,184:186])
dim(otu.notax) #1299  186
class(otu.notax) #dataframe
otu.notax <- as.matrix(otu.notax) #Turn 'otu.notax' into a matrix class
class(otu.notax) #matrix

#Create OTU table phyloseq object
OTU = otu_table(otu.notax, taxa_are_rows = TRUE)
head(OTU)
dim(OTU) #1299  186
class(OTU)
dim(otu.tax) #1299  188
head(otu.tax[,186:188])
tax.levels <- separate(data = otu.tax, 
                        col = Taxonomy, 
                        into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax.levels) #Notice that "Species" column is blank
dim(tax.levels) #1299  194
head(tax.levels[,186:194])
head(tax.levels[,188:193])
tax.only <- tax.levels[,188:193] #Keep only taxonomy columns from "Kingdom" up to "Genus"
head(tax.only)
dim(tax.only) #1299 6
class(tax.only) #data.frame
tax.m <- as.matrix(tax.only)
class(tax.m) #matrix
head(tax.m)

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax.m)
head(TAX)
dim(TAX) #1299 6
head(TAX)
class(TAX)

#Create phyloseq object 'phyloseqFlu' containing taxonomy, metadata, and OTU table
phyloseqFlu <- phyloseq(OTU, SAM, TAX)
phyloseqFlu #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1299 taxa and 196 samples ]
#sample_data() Sample Data:       [ 196 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1299 taxa by 6 taxonomic ranks ]

save(phyloseqFlu, file="SRD129Flu.phyloseq.RData")

############# ADONIS TEST FOR NASAL no bad samples, flu and control only, no DNEG12, DNEG6 ##############

#Run adonis function to determine effect of time and treatment on structure of nasal microbiota
adonis.Flu <- as(sample_data(phyloseqFlu), "data.frame")
class(adonis.Flu) #data.frame
dist.Flu <- distance(phyloseqFlu, method="bray") #Distance calculation using Bray-Curtis
set.seed(1) #Use set.seed function when running simulations to ensure all results are reproducible
full.Flu <- adonis(dist.Flu~Day*Treatment, data=adonis.Flu, permutations=9999)
full.Flu #Display results

#Call:
#adonis(formula = dist.Flu ~ Day * Treatment, data = adonis.Flu,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                   Df    SumsOfSqs MeanSqs   F.Model   R2        Pr(>F)    
#  Day              9     8.2853    0.92059   12.8444   0.37628   1e-04 ***
#  Treatment        1     0.3890    0.38903   5.4279    0.01767   1e-04 ***
#  Day:Treatment    9     1.4470    0.16077   2.2432    0.06571   1e-04 ***
#  Residuals        166   11.8977   0.07167             0.54034           
#  Total            185   22.0190                       1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#switched Treatment and Day order, obtained same conclusions: Day had the largest effect on variation
full.Flu_2 <- adonis(dist.Flu~Treatment*Day, data=adonis.Flu, permutations=9999)
full.Flu_2
#Output:
#Call:
#  adonis(formula = dist.Flu ~ Treatment * Day, data = adonis.Flu,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                   Df    SumsOfSqs MeanSqs   F.Model   R2        Pr(>F)    
#  Treatment        1     0.4068    0.40676   5.6752    0.01847   2e-04 ***
#  Day              9     8.2676    0.91862   12.8169   0.37548   1e-04 ***
#  Treatment:Day    9     1.4470    0.16077   2.2432    0.06571   1e-04 ***
#  Residuals        166   11.8977   0.07167             0.54034           
#  Total            185   22.0190                       1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#If distance function is giving an error message below:
#"Error: x should be a data.frame, data.table, tbl, tbl_df, array, or matrix."
#Use vegdist function from vegan package to run distance calculations instead of the distance function 
#(original "distance" function that was used below is no longer available) and use those calculations to run through adonis test
#vegdist requires that phyloseq object's OTU table has OTUs listed in the columns and sample names listed in rows. 
#Also, remove any OTUs with taxa_sums = 0 or non-numeric values. For example, this command can help remove OTUs with taxa_sums = 0: 
#OTU <- prune_taxa(taxa_sums(<yourOTUtable>) > 0, <yourOTUtable>)
#If you create a separate phyloseq object with this specific OTU table setup, 
#you should be able to run the vegdist function without any errors and use the output to run through adonis function

head(otu.meta4) #Sample names are listed in rows and OTUs are listed in columns in 'otu.meta4'
OTU.2 = otu_table(otu.meta4, taxa_are_rows = TRUE)
OTU.2.distance <- vegdist(OTU.2, method = "bray")
OTU.2.distance.full <- adonis(OTU.2.distance~Day*Treatment, data = adonis.Flu, permutations=9999)
OTU.2.distance.full
#Call:
#adonis(formula = OTU.2.distance ~ Day * Treatment, data = adonis.Flu,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Day             9    8.2853 0.92059 12.8444 0.37628  1e-04 ***
#  Treatment       1    0.3890 0.38903  5.4279 0.01767  1e-04 ***
#  Day:Treatment   9    1.4470 0.16077  2.2432 0.06571  1e-04 ***
#  Residuals     166   11.8977 0.07167         0.54034           
#Total         185   22.0190                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





#OLD CODE BELOW. IGNORE...


#############################-----Remove bad samples------############################################

#ATTEMPT #1 --- DID NOT WORK WHEN I TRIED RUNNING ADONIS FUNCTION
phyloseq.nw.nobad = subset_samples(phyloseq.nw, row.names(sample_data(phyloseq.nw)) != "SRD129_P150_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P151_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P152_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P153_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P154_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P156_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P163_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P164_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P165_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P166_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P167_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P168_N_D28")
#to remove specific samples, must include row.names(sample_data(phyloseq.nw)) != "<sample name>"
#add an "&" to add additional samples to exclude
phyloseq.nw.nobad
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 460 samples ]
#sample_data() Sample Data:       [ 460 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]
head(sample_data(phyloseq.nw.nobad))
tail(sample_data(phyloseq.nw.nobad))

#Save the phyloseq object without bad samples as a separate phyloseq object
save(phyloseq.nw.nobad, file="SRD129.phyloseq.nw.nobad.RData")


#ATTEMPT #2 --- ALSO DID NOT WORK ON ADONIS FUNCTION
#Read files
otunobad <- read.csv("SRD129abundsingleton2000OTUtablenobad.csv", row.names=1) #set column 1 as row names
metanobad <-read.csv("SRD129metadata06292018nobad.csv")

dim(otunobad) #2196 897
head(otunobad[,1:5])
head(otunobad[,890:897])

dim(metanobad) #892 7
head(metanobad[,1:7])

#Remove taxonomy
taxnobad <- otunobad[,(896:897)] #removed column 897 taxonomy and 896 to include the row names
head(taxnobad)
colnames(taxnobad)[1] <- "delete" #renamed column 1 (formerly 896) as "delete" which will later be deleted
head(taxnobad)

#OTU count data only
otunobad <- otunobad[,-897] #remove column 897 taxonomy to have only otu data
head(otunobad)
head(otunobad[,890:896]) 
dim(otunobad) #2196 rows 896 columns

#Transpose to match metadata format
otunobad.trans <- t(otunobad) #now rownames are sample names, columns are OTUs
head(otunobad.trans[,1:5])
head(metanobad) #row names are numbered, but want sample names as row names
metanobad <- metanobad[,-1]
dim(metanobad) #892 6
rownames(otunobad.trans) #contains 4 Mock samples

class(metanobad) #dataframe
class(otunobad) #dataframe

#Merge otu and meta data frames
otunobad.meta <- merge(metanobad, otunobad.trans, by.x=1, by.y=0) 
#x=1 means match via 1st column from meta; y=0 means match via rownames from otunobad.trans
head(otunobad.meta[,1:10])
tail(otunobad.meta)
class(otunobad.meta) #dataframe
dim(otunobad.meta) #892 2202

#With "All" column
write.csv(otunobad.meta, file="SRD129abundsingleton2000.otu.meta_allnobad.csv")

#Subset tissues with "All" column
nw.nobad <- subset(otunobad.meta, Tissue=="N")
head(nw.nobad[,1:10])
dim(nw.nobad) #460 rows 2202 columns

write.csv(nw.nobad, file="SRD129abundsingleton2000.otu.meta_allnobad.nasalonly.csv")

#########-------NASAL PHYLOSEQ with "All" column, no bad samples, flu and control only------##########
#########-------no DNEG12, DNEG6-------#########

#Nasal
#Split info again
#meta only
head(nw.nobad[,1:10])
dim(nw.nobad) #460 2202
colnames(nw.nobad)
rownames(nw.nobad)
nw.nobad<- nw.nobad[!grepl("PRRSV", nw.nobad$Treatment),]
nw.nobad<- nw.nobad[!grepl("BB", nw.nobad$Treatment),]
nw.nobad <- nw.nobad[!grepl("DNEG12", nw.nobad$Day),]
nw.nobad <- nw.nobad[!grepl("DNEG6", nw.nobad$Day),]
unique(nw.nobad$Day) #D0  D1  D10 D14 D21 D3  D36 D42 D7  D28
unique(nw.nobad$Treatment) #IAV Control

meta.nw.nobad <- nw.nobad[,1:6] #take columns 1-6 to make meta.nw.nobad
head(meta.nw.nobad)
row.names(meta.nw.nobad) <- meta.nw.nobad[,1] #make column 1 be rownames for meta.nw.nobad
head(meta.nw.nobad)
meta.nw.nobad <- meta.nw.nobad[,-1] #remove the extra "Sample" column
head(meta.nw.nobad)
dim(meta.nw.nobad) #186 5
class(meta.nw.nobad) #data.frame


#Create SAM metada table phyloseq object
SAMnw.nobad = sample_data(meta.nw.nobad, errorIfNULL = TRUE)
head(SAMnw.nobad)
dim(SAMnw.nobad) #186 5

#OTU only
head(nw.nobad[,1:10])
dim(nw.nobad) #186 2202
head(nw.nobad[,2199:2202])
otu.nw.nobad <- nw.nobad[,c(1,7:2202)] #select Sample and Otu columns to create otu.nw.nobad dataframe
head(otu.nw.nobad[,1:10])
dim(otu.nw.nobad) #186 2197
row.names(otu.nw.nobad) <- otu.nw.nobad[,1] #make first column be rownames for otu.nw.nobad
head(otu.nw.nobad[,1:10])
otu.nw.nobad <- otu.nw.nobad[,-1] #remove the first column
head(otu.nw.nobad[,1:10])
otu.nw.nobad.trans <- t(otu.nw.nobad) #transpose otu.nw.nobad to have Otu as rownames, sample names as column names
head(otu.nw.nobad.trans[,1:10])
dim(otu.nw.nobad.trans) #2196 186
head(otu.nw.nobad.trans[,180:186])

#Merge tax back into otu for correct format and taxons
head(taxnobad)
otu.tax.nw.nobad <- merge(otu.nw.nobad.trans, tax, by=0) #merge by rownames aka Otu rownames
dim(otu.tax.nw.nobad) #2196 2189
head(otu.tax.nw.nobad[,185:189])
head(otu.tax.nw.nobad[,1:5])
row.names(otu.tax.nw.nobad) <- otu.tax.nw.nobad[,1] #set first row as rownames
head(otu.tax.nw.nobad[,1:5])
otu.tax.nw.nobad <- otu.tax.nw.nobad[,-1] #remove first row, extraneous Otu column
head(otu.tax.nw.nobad[,1:5])
dim(otu.tax.nw.nobad) #2196 188
head(otu.tax.nw.nobad[,185:188])

#Split again
dim(otu.tax.nw.nobad) #2196 188
head(otu.tax.nw.nobad[,185:188])
otu.notax.nw.nobad <- otu.tax.nw.nobad[,1:186] #take rows 1-186 to make new dataframe otu.notax.nw (187 is delete column, 188 is taxonomy column)
head(otu.notax.nw.nobad[,1:5])
head(otu.notax.nw.nobad[,180:186])
dim(otu.notax.nw.nobad) #2196 186
class(otu.notax.nw.nobad) #dataframe
otu.notax.nw.nobad <- as.matrix(otu.notax.nw.nobad) #turn otu.notax.nw.nobad into a matrix class
class(otu.notax.nw.nobad) #matrix

#Create OTU table phyloseq object
OTUnw.nobad = otu_table(otu.notax.nw.nobad, taxa_are_rows=TRUE)
head(OTUnw.nobad)
dim(OTUnw.nobad) #2196 186
class(OTUnw.nobad)

dim(otu.tax.nw.nobad) #2196 188
head(otu.tax.nw.nobad[,185:188])
tax.levels.nw.nobad <- separate(data = otu.tax.nw.nobad, 
                                col = Taxonomy, 
                                into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled Kingdom to Species
head(tax.levels.nw.nobad) #notice that Species column is blank
dim(tax.levels.nw.nobad) #2196 194
head(tax.levels.nw.nobad[,187:194])
tax.only.nw.nobad <- tax.levels.nw.nobad[,188:193] #keep only taxonomy columns Kingdom to Genus
head(tax.only.nw.nobad)
dim(tax.only.nw.nobad) #2196 6
class(tax.only.nw.nobad) #data.frame
tax.m.nw.nobad <- as.matrix(tax.only.nw.nobad)
class(tax.m.nw.nobad) #matrix
head(tax.m.nw.nobad)
#tax.m.nw.nobad <- tax.m.nw.nobad[,-1]

#Create TAX taxonomy table phyloseq object
TAXnw.nobad = tax_table(tax.m.nw.nobad)
head(TAXnw.nobad)
dim(TAXnw.nobad) #2196 6
class(TAXnw.nobad)

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.nw.nobad <- phyloseq(OTUnw.nobad, SAMnw.nobad, TAXnw.nobad)
phyloseq.nw.nobad #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 186 samples ]
#sample_data() Sample Data:       [ 186 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]

save(phyloseq.nw.nobad, file="SRD129.phyloseq.nw.nobad2.RData")

##############################################################################################################################

############# ADONIS TEST FOR NASAL no bad samples, flu and control only, no DNEG12, DNEG6 ##############

adonis.nw.nobad <- as(sample_data(phyloseq.nw.nobad), "data.frame")
class(adonis.nw.nobad) #data.frame
dist.nw.nobad <- distance(phyloseq.nw.nobad, method="bray")
set.seed(1)
full.nw.nobad <- adonis(dist.nw.nobad~Day*Treatment, data=adonis.nw.nobad, permutations=9999)
full.nw.nobad

#Call:
#adonis(formula = dist.nw.nobad ~ Day * Treatment, data = adonis.nw.nobad,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Day             9    8.2429 0.91588 12.7536 0.37484  1e-04 ***
#Treatment       1    0.3869 0.38688  5.3872 0.01759  1e-04 ***
#Day:Treatment   9    1.4394 0.15994  2.2271 0.06546  2e-04 ***
#Residuals     166   11.9211 0.07181         0.54211           
#Total         185   21.9903                 1.00000 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1