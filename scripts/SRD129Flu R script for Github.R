#R scripts for analyzing and visualizing data from the research paper 
"Changes in the swine nasal microbiota following influenza A virus challenge in a longitudinal study"

#################################################################################################################################################################################################

#SECOND SECTION (after processing sequences in mothur): Creating OTU Table

#Purpose: Create OTU table with R using specific files generated from mothur

#Files needed:
#SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#SRD129Flu.outsingletons.abund.taxonomy.csv
#SRD129Flu.outsingletons.abund.opti_mcc.0.03.subsample.shared
#SRD129Flu.outsingletons.abund.subsample.shared.csv

#The output files from mothur needed in R for this section and subsequent sections, aside from fasta files, are text files that can be saved as csv for ease of use in R.

#Load library package
library(tidyverse)

#To start creating OTU table, edit taxonomy file
#Save "SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy" file generated from mothur as a csv file in a spreadsheet editor. 
#Named file as "SRD129Flu.outsingletons.abund.taxonomy.csv"
taxonomy <- read.csv("./data/SRD129Flu.outsingletons.abund.taxonomy.csv") #Import this csv file from working directory using "read.csv" function
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Substitute (100) and variations of that with nothing ('') from the Taxonomy column
taxonomy[1:6,3] #Show column number 3, rows 1 through 6 in 'taxonomy' dataframe
write.csv(taxonomy, file = "SRD129abundsingleton2000taxonomy.csv") #Write 'taxonomy' into a csv file and open the csv file in a spreadsheet editor to remove the size and numbered rows

#Edit subsample.shared file
#Save "SRD129Flu.outsingletons.abund.opti_mcc.0.03.subsample.shared" file generated from mothur as a csv file in a spreadsheet editor. 
#Named file as "SRD129Flu.outsingletons.abund.subsample.shared.csv"
shared <- read.csv("./data/SRD129Flu.outsingletons.abund.subsample.shared.csv", stringsAsFactors = FALSE)
head(shared) #Check on the first part of 'shared' dataframe
shared[1:6,1]
shared <- t(shared) #Transpose 'shared'
head(shared)
shared[1:6,1]
write.csv(shared, file = 'SRD129abundsingleton2000shared.csv') #Open this csv file in a spreadsheet editor and remove the "V*", "label", and "numOtus" rows

#Read the revised SRD129abundsingleton2000shared.csv file
shared <- read.csv("./data/SRD129abundsingleton2000shared.csv")
colnames(shared) [1] <- "OTU" #Rename first column of 'shared' to "OTU"
taxonomy <- read.csv("./data/SRD129abundsingleton2000taxonomy.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge 'shared' and 'taxonomy' objects by OTU
head(OTUtable)
nrow(OTUtable) #Count number of rows in 'OTUtable'
ncol(OTUtable) #Count number of columns in 'OTUtable'
write.csv(OTUtable, file= "SRD129abundsingleton2000OTUtable.csv") 
#Open this csv file in a spreadsheet editor and remove the first row (numbered rows) and "Size" column

#Check OTU table
OTUtable <-read.csv("./data/SRD129abundsingleton2000OTUtable.csv", stringsAsFactors = FALSE)
head(OTUtable) #Check the first part of 'OTUtable' to make sure it looks ok

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

#################################################################################################################################################################################################

#THIRD SECTION: Creating phyloseq objects for each tissue

#Purpose: Create phyloseq objects to be used to calculate alpha and beta diversity measures for nasal and tonsil tissue samples. This section will also use the adonis function to determine the effect of time and treatment on the community structure of nasal and tonsil microbiota.

#Files needed:
#SRD129abundsingleton2000OTUtable.csv
#SRD129metadata.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(vegdist)

#Read files
otu <- read.csv("./data/SRD129abundsingleton2000OTUtable.csv", row.names=1) #Set column 1 as row names
meta <-read.csv("./data/SRD129metadata.csv")

dim(otu) #Check dimensions of 'otu'
head(otu[,1:5]) #Check the first part of 'otu' table
head(otu[,185:187]) #Check last part of 'otu' table
dim(meta) #Determine the dimensions of 'meta' dataframe. It has 186 rows and 5 columns
head(meta[,1:5]) #Check first part of 'meta' table

#Remove taxonomy from 'otu'
tax <- otu[,(186:187)] #Remove column 186 taxonomy and 187 to copy the row names from 'otu' to 'tax'
head(tax)
colnames(tax)[1] <- "delete" #Rename column 1 of 'tax' (formerly column 186) as "delete" which will be deleted later
head(tax)

#Modify 'otu' with only OTU count data
otu <- otu[,-187] #Remove column 187 taxonomy in 'otu' to have only OTU data
head(otu[,180:186]) 
dim(otu) #Dimensions of 'otu' show 1299 rows 186 columns

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu) 
#Now rownames in 'otu.trans' are sample names, columns are OTUs
head(otu.trans[,1:5])
head(meta) #Row names are numbered, but we want sample names as row names
rownames(meta) <- meta$Sample #Set names in "Sample" as rownames in 'meta'
head(meta)
class(meta) #The type of class that 'meta' is is dataframe
class(otu) #dataframe

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) #Merge by names of the columns that are common to both x and y (columns with common names between the two data sets)
#by.x=0 means match by rownames in 'meta'; y=0 means match by rownames in 'otu.trans'
head(otu.meta[,1:10])
class(otu.meta) #Check class type of 'otu.meta'. It should be a dataframe.

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2' dataframe
otu.meta2<- cbind(otu.meta) #Make second copy of 'otu.meta' to use to include an "All" column
otu.meta2$All <- with(otu.meta2, paste0(Day, sep=" ", Treatment)) #Combine "Day" and "Treatment" columns into an "All" column
head(otu.meta2) #Check first part of otu.meta2
dim(otu.meta2) #Check dimensions of 'otu.meta2' dataframe
head(otu.meta2[,1300:1305]) #Check the first part of the end of 'otu.meta2' dataframe
rownames(otu.meta2) <- otu.meta2$Row.names #Set "Row.names" as rownames
otu.meta2 <- otu.meta2[,-1]
head(otu.meta2[,1:10])
dim(otu.meta2) #186 1304
head(otu.meta2[,1300:1304])
otu.meta2<- otu.meta2[,c(1:4,1304,5:1303)] #Reorder columns to have "All" column after "Treatment" column
head(otu.meta2[,1:10]) #Check the first part of the beginning of 'otu.meta2' dataframe
head(otu.meta2[,1300:1303])
write.csv(otu.meta2, file="SRD129abundsingleton2000.otu.meta.csv")

#Creating phyloseq objects for nasal samples

#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
dim(otu.meta2) #186 1643
otu.meta3 <- otu.meta2[,1:5] #Take columns 1-5 to make 'otu.meta3'
head(otu.meta3)
dim(otu.meta3) #186 5

#Create SAM metadata table phyloseq object
SAM = sample_data(otu.meta3, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #186 5

#Pull out OTU data from 'otu.meta2' dataframe
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

#Split 'otu.tax.nw2' again
dim(otu.tax) #1299 188
head(otu.tax[185:188])
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


#Run adonis function to determine effect of time and treatment on structure of nasal microbiota

adonis.Flu <- as(sample_data(phyloseqFlu), "data.frame")
class(adonis.Flu) #data.frame

#distance function to run distance calculations
dist.Flu <- distance(phyloseqFlu, method="bray") #Distance calculation using Bray-Curtis
set.seed(1) #Use set.seed function when running simulations to ensure all results are reproducible
full.Flu <- adonis(dist.Flu~Day*Treatment, data=adonis.Flu, permutations=9999)
full.Flu #Display results
#Output:
#Call:
#Call:
#  adonis(formula = dist.Flu ~ Day * Treatment, data = adonis.Flu,      permutations = 9999) 

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
    
###Day had largest effect on variation
    
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
#Output is the same as full.Flu:
#Call:
#adonis(formula = OTU.2.distance ~ Day * Treatment, data = adonis.Flu,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                   Df    SumsOfSqs MeanSqs   F.Model   R2        Pr(>F)    
#  Day              9     8.2853    0.92059   12.8444   0.37628   1e-04 ***
#  Treatment        1     0.3890    0.38903   5.4279    0.01767   1e-04 ***
#  Day:Treatment    9     1.4470    0.16077   2.2432    0.06571   1e-04 ***
#  Residuals      166     11.8977   0.07167             0.54034           
#  Total          185     22.0190                       1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#################################################################################################################################################################################################

#FOURTH SECTION: Beta diversity
    
#Purpose: This code generates non-metric multidimensional scaling ordination based on Bray-Curtis dissimilarities to create NMDS plots, and runs pairwise.adonis function to identify any significant differences in bacterial composition between treatment groups on a given day.
    
#Files needed:
#phyloseqFlu

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
    
#Load this function from the funfuns R package (https://github.com/Jtrachsel/funfuns)
NMDS_ellipse <- function(metadata, OTU_table, grouping_set,
                         distance_method = 'bray',
                         rand_seed = 77777,
                         MDS_trymax = 1000,
                         autotransform = FALSE,
                         wascores = TRUE,
                         expand = FALSE){
  require(vegan)
  require(tidyr)
  
  if (grouping_set %in% colnames(metadata)){
    if (all(rownames(metadata) == rownames(OTU_table))){
      
      set.seed(rand_seed)
      generic_MDS <- metaMDS(OTU_table, k = 2,
                             trymax = MDS_trymax,
                             autotransform = autotransform,
                             distance = distance_method,
                             wascores = wascores,
                             expand = expand)
      
      stress <- generic_MDS$stress
      nmds_points <- as.data.frame(generic_MDS$points)
      metadata <- metadata[match(rownames(generic_MDS$points), rownames(metadata)),]
      metadata <- as.data.frame(metadata) # weird things were happening when a grouped tibble was used as metadata...
      metanmds <- cbind(metadata, nmds_points)
      # browser()
      nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
      metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function
      
      #check to make sure at least 3 obs for each grouping_set
      
      numobs <- metanmds %>% group_by(!!grouping_set) %>% summarise(n=n())
      if (min(numobs$n) >= 3){
        ord <- ordiellipse(generic_MDS, metanmds[[grouping_set]], label = TRUE, conf = .95, kind = 'se', draw = 'none')
        
        df_ell <- data.frame()
        for (d in levels(metanmds[[grouping_set]])){
          df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds[[grouping_set]] == d,],
                                                           vegan:::veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
        }
        
        
        
        # this loop assigns the group centroid X coordinates to each sample, there is probably a better way...
        
        metanmds$centroidX <- NA
        metanmds$centroidY <- NA
        
        
        for (level in levels(metanmds[[grouping_set]])){
          metanmds[metanmds[[grouping_set]] == level,]$centroidX <- nmds.mean$MDS1[nmds.mean$group == level]
          metanmds[metanmds[[grouping_set]] == level,]$centroidY <- nmds.mean$MDS2[nmds.mean$group == level]
          
          
        }
        print(paste('Ordination stress:', stress, sep = ' '))
        return(list(metanmds, df_ell, generic_MDS))
        
      } else {
        warning('One of your groups in "grouping_set" has less than 3 observations, cannot generate elipses')
        df_ell <- data.frame()
        return(list(metanmds, df_ell, generic_MDS))}
      
    } else {
      stop('The rownames for your OTU table and metadata do not match.')
    }
    
  }else {
    stop('The grouping set column you have provided in not in your metadata.')
  }
  
  
  
}


#Load the following pairwise.adonis function, taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

    
#Setting up 'phyloseqFlu' into dataframes for NMDS calculation
flu.sam <- data.frame(phyloseqFlu@sam_data) #Make 'phyloseqFlu sam_data' into dataframe
flu.otu <- data.frame(t(phyloseqFlu@otu_table)) #Make 'phyloseqFlu otu_table' into dataframe
class(flu.sam) #data.frame
rownames(flu.sam) == row.names(flu.otu) #For rows with sums greater than 1 in 'flu.otu', move rows and their respective sum values into "numOTUs" column in 'flu.sam'
flu.sam$numOTUS <- rowSums(flu.otu > 1)
head(flu.sam)

#NMDS calculation
flu.otu[1:10,1:10]
flu.NMDS <- NMDS_ellipse(flu.sam, flu.otu, grouping_set = 'All')
#Output:
#Result: [1] "Ordination stress: 0.168876156130503"
    
#Separate meta data and ellipse data to two lists to make NMDS plot
flu.NMDS.2 <- flu.NMDS[[1]] 
#'flu.NMDS.2' has meta data + MDS calculations. Select this 1st list of 'flu.NMDS' using double brackets
flu.df_ell.2 <- flu.NMDS[[2]]
#'flu.df_ell.2' is accessing 2nd list from 'flu.NMDS' that has ellipse calculations
#Need two separate lists for making NMDS plot
flu.df_ell.2$group #20 levels
head(flu.df_ell.2)

#Create "Day" and "Treatment" columns within 'nw.df_ell' for faceting purposes
flu.df_ell.2$Day <- sub(' [A-Za-z]+', '\\1', flu.df_ell.2$group) #Created "Day" column, '\\1' returns the first part of the regular expression D[A-z0-9]+ from 'flu.df_ell.2$group'
flu.df_ell.2$Treatment <- sub('D[A-z0-9]+ ', '\\2', flu.df_ell.2$group) #Create "Treatment" column, '\\2' returns the second part of the sub expression ([A-Za-z]+) from 'flu.df_ell.2$group'
head(flu.df_ell.2)

#Restructure level order for 'nw.metanmds' and 'nw.df_ell'
flu.NMDS.2$Day = factor(flu.NMDS.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
flu.df_ell.2$Day = factor(flu.df_ell.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(flu.df_ell.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
levels(flu.NMDS.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
    
#Creating plot from NMDS calculations
FluNMDSplot <- ggplot(data=flu.NMDS.2, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
  geom_path(data=flu.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15)) +
  labs(caption = 'Ordination stress = 0.169', color="Treatment group")
FluNMDSplot
    
#Save 'FluNMDSplot' as a .tiff for publication, 500dpi
ggsave("Figure_1.tiff", plot=FluNMDSplot, width = 10, height = 6, dpi = 500, units =c("in"))
    
#Using pairwise.adonis function
flu.adon <- pairwise.adonis(flu.otu, flu.sam$All, sim.method = 'bray', p.adjust.m = 'bonferroni') 
#Run pairwise.adonis on 'flu.otu' OTU table and "All" column of 'flu.sam' dataframe
flu.adon$pairs #List all comparisons in the "pairs" column of 'flu.adon'
goodcomps.flu.adon <- c(grep('D0 [A-Za-z]+ vs D0 [A-Za-z]+', flu.adon$pairs),
                        grep('D1 [A-Za-z]+ vs D1 [A-Za-z]+', flu.adon$pairs),
                        grep('D3 [A-Za-z]+ vs D3 [A-Za-z]+', flu.adon$pairs),
                        grep('D7 [A-Za-z]+ vs D7 [A-Za-z]+', flu.adon$pairs),
                        grep('D10 [A-Za-z]+ vs D10 [A-Za-z]+', flu.adon$pairs),
                        grep('D14 [A-Za-z]+ vs D14 [A-Za-z]+', flu.adon$pairs),
                        grep('D21 [A-Za-z]+ vs D21 [A-Za-z]+', flu.adon$pairs),
                        grep('D28 [A-Za-z]+ vs D28 [A-Za-z]+', flu.adon$pairs),
                        grep('D36 [A-Za-z]+ vs D36 [A-Za-z]+', flu.adon$pairs),
                        grep('D42 [A-Za-z]+ vs D42 [A-Za-z]+', flu.adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
goodcomps.flu.adon
goodcomps.flu.adon.2 <- flu.adon[goodcomps.flu.adon,] #Rename 'goodcomps.flu.adon' vector to 'goodcomps.flu.adon.2'
goodcomps.flu.adon.2
goodcomps.flu.adon.2$p.adjusted <- p.adjust(goodcomps.flu.adon.2$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
goodcomps.flu.adon.2$p.adjusted2 <- round(goodcomps.flu.adon.2$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
goodcomps.flu.adon.2$p.adjusted2[goodcomps.flu.adon.2$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(goodcomps.flu.adon.2, file='Flu.pairwisecomparisons.csv', row.names=TRUE)

#################################################################################################################################################################################################

#FIFTH SECTION: Alpha diversity
    
#Purpose: This code calculates alpha diversity metrics (Shannon, Inverse Simpson) that are plotted as box and whisker plots, and uses wilcoxon rank sum test to assess any significant differences in diversity between treatment groups on a given day.

#Files needed:
#phyloseqFlu

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
    
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
flu.sam$shannon <- diversity(otu.meta4) #"diversity" is a vegan function. The default index is set at "shannon". I added a shannon index column in 'flu.sam'
flu.sam$invsimpson <- diversity(otu.meta4,index = 'invsimpson') #We used 'invsimpson' since it is easier to interpret than Simpson values and won't need to "inverse" the Simpson values to understand (With Simpson values, the lower the number, the higher the diversity)
flu.sam$Day = factor(flu.sam$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42")) # Set the level order of values in "Day" column
levels(sample_data(flu.sam)$Day) #Set the level order of values in "Day" column to D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"   
    
#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within flu.sam
average.shannon.invsimpson.numOTUs <- aggregate(flu.sam[, 6:8], list(flu.sam$All), mean)
print(average.shannon.invsimpson.numOTUs)
#Output:
#       Group.1   numOTUS  shannon invsimpson
#1   D0 Control  69.10000 2.565335   6.377505
#2       D0 IAV  51.20000 2.169701   4.399913
#3   D1 Control  50.90000 2.253141   5.126906
#4       D1 IAV  45.20000 2.030144   4.128696
#5  D10 Control  62.70000 2.445356   5.112930
#6      D10 IAV  51.70000 2.305930   5.136920
#7  D14 Control  81.10000 2.853185   7.046119
#8      D14 IAV  92.20000 3.069216   8.757443
#9  D21 Control  79.90000 2.925234   7.111600
#10     D21 IAV  91.33333 3.135175   8.865022
#11 D28 Control 105.50000 3.511063  13.981237
#12     D28 IAV 103.20000 3.608026  15.135310
#13  D3 Control  85.60000 2.902182   7.333515
#14      D3 IAV  68.00000 2.571588   5.847328
#15 D36 Control  91.80000 3.152522  13.080790
#16     D36 IAV 109.11111 3.409859  15.417720
#17 D42 Control 108.40000 3.515769  17.795337
#18     D42 IAV 113.77778 3.477660  16.745868
#19  D7 Control  56.80000 2.382253   5.533712
#20      D7 IAV  26.20000 1.817367   4.238963
write.csv(average.shannon.invsimpson.numOTUs, file="Flu.average.shannon.invsimpson.num.OTUs.txt", row.names=TRUE)
flu.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(flu.sam$shannon, flu.sam$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(flu.pairwise.wilcox.shannon.test) #Look at the results of 'Flu.pairwise.wilcox.shannon.test'
flu.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(flu.sam$invsimpson, flu.sam$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the inverse simpson indices in "invsimpson" column
print(flu.pairwise.wilcox.invsimpson.test)
    
#Generate a box and whisker plot of shannon
flu.shannon <- ggplot(data = flu.sam, aes(x=Treatment, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
# "free" within "facet_wrap" allows each plot to customize the scale to the specific data set (no forced scaling applied to all plots)
# "position = position_dodge2(preserve = 'total')" fixes the ggplot box width, making them wider, prevents narrow boxes from forming in the plot
flu.shannon
    
#Generate a box and whisker plot of inverse simpson diversity
flu.invsimp <- ggplot(data = flu.sam, aes(x=Treatment, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
flu.invsimp
    
#Combining plots
flu.combalpha <- plot_grid(flu.shannon + theme(legend.position = "none"), flu.invsimp, labels = "AUTO")
flu.combalpha

#Save 'flu.combalpha' as a .tiff for publication, 500dpi
ggsave("Figure_4.tiff", plot=flu.combalpha, width = 10, height = 5, dpi = 500, units =c("in"))

#################################################################################################################################################################################################

#SIXTH SECTION: Magnitude of Change in Nasal Microbiota
    
#Purpose: This code plots the F-statistic from PERMANOVA pairwise comparisons of Control and IAV groups, over time (displays the magnitude of change in the nasal bacterial community structure of the IAV group relative to control)
    
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
    
#The F-values were obtained from 'Flu.pairwisecomparisons.csv' data generated in the "Beta diversity" section
#Flu_MagnitudeOfChange.csv was created from Flu.pairwisecomparisons.csv
#To make the Flu_MagnitudeOfChange.csv file, open Flu.pairwisecomparisons.csv file in excel, copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. Add "Day" and "Treatment" columns and save as "Flu_MagnitudeOfChange.csv".
    
flu.mag <- read.csv("./data/Flu_MagnitudeOfChange.csv")
class(flu.mag)
flu.mag$Day <- factor(flu.mag$Day) #Encode "Day" as a factor
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
    
#################################################################################################################################################################################################

#SEVENTH SECTION: Differential Abundance of Genera in Nasal Microbiota using DESeq2
    
#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially abundant between treatment groups and control group
    
#Files needed:
#Mothur shared file: FS1bfinal.outsingletons.abund.opti_mcc.shared
#Mothur constaxonomy file: FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#Metadata: FS1babundsingleton2000metadata.csv
#FS1b_FinalDiffAbundNasalGenus_IC.csv
#FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library(cowplot)
    
#Annotations
#fc = infeed, control
#ic = injected, control
#fi = infeed, injected
#if = injected, infeed
    
#Preparing objects for DESeq2: load files
otu2 <- import_mothur(mothur_shared_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.shared')
taxo2 <- import_mothur(mothur_constaxonomy_file = './data/SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
taxo2
meta2 <- read.table(file = './data/SRD129metadata.csv', sep = ',', header = TRUE)
    
#Organize 'meta2' meta file
rownames(meta2) <- meta2$Sample #Make names in "Sample" become rownames for 'meta2'
meta2 <- meta2[,-1] #Remove Sample column
meta2$Set <- paste(meta2$Day, meta2$Treatment, sep = '_')

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta2 <- sample_data(meta2) 
SRD129 <- phyloseq(otu2, taxo2)
SRD129 <- merge_phyloseq(SRD129, phy_meta2) #Combines the 'phy_meta2' metadata with this phyloseq object
colnames(tax_table(SRD129))
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
SRD129


####Continue here##
    
#Prune samples from 'FS1b'
FS1b <- prune_samples(sample_sums(FS1b) > 2000, FS1b)  #This removes samples that have fewer than 2000 sequences associated with them.
FS1b <- prune_taxa(taxa_sums(FS1b) > 10, FS1b)        #This removes OTUs that occur less than 10 times globally
tax_table(FS1b) [1:5, 1:6] #Checking what's in 'tax_table' first 5 rows, first 6 columns
    
#Grouping OTUs by Genus using the tax_glom function 
FS1b.genus <- tax_glom(FS1b, taxrank = "Genus")
#This method merges species that have the same taxonomy at a certain taxanomic rank (in this case, by "Genus"). 
#Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 
    
######################################## PRIMARY COMPARISONS TO MAKE ####################################################
    
# NMDS plot showed that disperion is different between days, so I subsetted by day and tissue
    
# Important comparisons to make (significant changes in beta diversity between treatments)
    
# Compare Days 4, 7, 11 NON compared to each of the two treatments 
# Compare Days 4, 7, 11 IM and IF
# Compare Day 14 NON and IF
    
# Other comparisons to make (no significant changes in beta diversity between treatments)
    
# Compare Days 0, 14 IM and IF
# Compare Day 0 NON compared to each of the two treatments
# Compare Day 14 NON and IM
    
##################################################### Day 0 ############################################################
    
############## Day 0 Nasal #########################
    
sample_data(FS1b.genus)
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D0.nw.De'
FS1b.D0.nw <- subset_samples(FS1b.genus, Day == 0 & Tissue == 'NW')
sample_sums(FS1b.D0.nw) #Returns the total number of individuals observed from each sample
FS1b.D0.nw <- prune_taxa(taxa_sums(FS1b.D0.nw) > 1, FS1b.D0.nw) #If taxa_sums is >1, then it will print that out in 'FS1b.D0.nw' object and not include any samples with taxa of sums <1.
rowSums(FS1b.D0.nw@otu_table)
FS1b.D0.nw.De <- phyloseq_to_deseq2(FS1b.D0.nw, ~ Set) # "~Set" : define the major sample covariate as the study design factor. This will be whatever you want to group data by, whatever column you used to designate ellipses with
FS1b.D0.nw.De <- DESeq(FS1b.D0.nw.De, test = "Wald", fitType = "parametric") #Differential expression analysis based on the negative binomial (aka Gamma-Poisson) distribution
    
######### Day 0 Nasal IF vs NON ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "0_NW_IF")
#Output:
#IF = 19
sum(meta2$Set == "0_NW_NON")
#Output:
#NON = 18
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.fc' table
res.D0.fc = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IF","0_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D0.fc
sigtab.D0.fc = res.D0.fc[which(res.D0.fc$padj < .05), ]
sigtab.D0.fc = cbind(as(sigtab.D0.fc, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.fc), ], "matrix")) #"cbind" combines all columns together, regardless of rownames (if you want to match by rownames, use merge function)
format(sigtab.D0.fc$padj, scientific = TRUE)
sigtab.D0.fc$newp <- format(round(sigtab.D0.fc$padj, digits = 3), scientific = TRUE)
sigtab.D0.fc$Treatment <- ifelse(sigtab.D0.fc$log2FoldChange >=0, "IF", "NON") #Assigning "IF" = yes, "NON" = no; it's important to make sure you have the correct group names in the "yes" and "no" position for "ifelse" function
sigtab.D0.fc
    
#Summarize 'sigtab.D0.fc' table
sum.sigtab.D0.fc <- summary(sigtab.D0.fc)
sum.sigtab.D0.fc
    
#plot 'sigtab.D0.fc'
deseq.D0.fc <- ggplot(sigtab.D0.fc, aes(x=reorder(rownames(sigtab.D0.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 0')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
# "reorder" makes the logfold changes appear in numerical order (largest logfold value at the top and at the bottom of the plot), making it easier to see
deseq.D0.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.fc'
sigtab.D0.fc$OTU <- rownames(sigtab.D0.fc)
sigtab.D0.fc$comp <- 'D0_nasal_IFvsNON'
    
#Create a final table ('final.nonsigtab') that lists all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D0.fc'. 
#Other within-day comparisons that had no significant changes in beta diversity between two treatment groups will be added to 'final.nonsigtab'.
final.nonsigtab <- sigtab.D0.fc
    
########## Day 0 Nasal IM vs NON  ####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "0_NW_IM")
#Output:
#IM = 19
sum(meta2$Set == "0_NW_NON")
#Output:
#NON = 18
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.ic' table
FS1b.D0.nw.De$Set
res.D0.ic = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IM","0_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0.ic = res.D0.ic[which(res.D0.ic$padj < .05), ]
sigtab.D0.ic = cbind(as(sigtab.D0.ic, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.ic), ], "matrix"))
format(sigtab.D0.ic$padj, scientific = TRUE)
sigtab.D0.ic$newp <- format(round(sigtab.D0.ic$padj, digits = 3), scientific = TRUE)
sigtab.D0.ic$Treatment <- ifelse(sigtab.D0.ic$log2FoldChange >=0, "IM", "NON")
    
#Summarize 'sigtab.D0.ic' table
sum.sigtab.D0.ic <- summary(sigtab.D0.ic)
sum.sigtab.D0.ic
    
#Plot 'sigtab.D0.ic'
deseq.D0.ic <- ggplot(sigtab.D0.ic, aes(x=reorder(rownames(sigtab.D0.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 0')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4E9', '#999999'))
deseq.D0.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.ic'
sigtab.D0.ic$OTU <- rownames(sigtab.D0.ic)
sigtab.D0.ic
sigtab.D0.ic$comp <- 'D0_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D0.ic' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.ic)
    
######### Day 0 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "0_NW_IF")
#Output:
#IF = 19
sum(meta2$Set == "0_NW_IM")
#Output:
#IM = 19
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.if' table
res.D0.if = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IM","0_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D0.if
sigtab.D0.if = res.D0.if[which(res.D0.if$padj < .05), ]
sigtab.D0.if = cbind(as(sigtab.D0.if, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.if), ], "matrix"))
format(sigtab.D0.if$padj, scientific = TRUE)
sigtab.D0.if$newp <- format(round(sigtab.D0.if$padj, digits = 3), scientific = TRUE)
sigtab.D0.if$Treatment <- ifelse(sigtab.D0.if$log2FoldChange >=0, "IM", "IF")
sigtab.D0.if
    
#Summarize 'sigtab.D0.if' table
sum.sigtab.D0.if <- summary(sigtab.D0.if)
sum.sigtab.D0.if
    
#Plot 'sigtab.D0.if'
deseq.D0.if <- ggplot(sigtab.D0.if, aes(x=reorder(rownames(sigtab.D0.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 0')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust = 0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "IM"), values = c('#E69F00', '#56B4E9'))
deseq.D0.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.if'
sigtab.D0.if$OTU <- rownames(sigtab.D0.if)
sigtab.D0.if$comp <- 'D0_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D0.if' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.if)
    
    
    
##################################################### Day 4 ############################################################
    
############## Day 4 Nasal #########################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D4.nw.De'
FS1b.D4.nw <- subset_samples(FS1b.genus, Day == 4 & Tissue == 'NW')
sample_sums(FS1b.D4.nw)
FS1b.D4.nw <- prune_taxa(taxa_sums(FS1b.D4.nw) > 1, FS1b.D4.nw)
rowSums(FS1b.D4.nw@otu_table)
FS1b.D4.nw.De <- phyloseq_to_deseq2(FS1b.D4.nw, ~ Set)
FS1b.D4.nw.De <- DESeq(FS1b.D4.nw.De, test = "Wald", fitType = "parametric")
    
######### Day 4 Nasal IF vs NON ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "4_NW_IF")
#Output:
#IF = 20
sum(meta2$Set == "4_NW_NON")
#Output:
#NON = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.fc' table
res.D4.fc = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IF","4_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D4.fc
sigtab.D4.fc = res.D4.fc[which(res.D4.fc$padj < .05), ]
sigtab.D4.fc = cbind(as(sigtab.D4.fc, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.fc), ], "matrix"))
format(sigtab.D4.fc$padj, scientific = TRUE)
sigtab.D4.fc$newp <- format(round(sigtab.D4.fc$padj, digits = 3), scientific = TRUE)
sigtab.D4.fc$Treatment <- ifelse(sigtab.D4.fc$log2FoldChange >=0, "IF", "NON")
sigtab.D4.fc
    
#Plot 'sigtab.D4.fc'
deseq.D4.fc <- ggplot(sigtab.D4.fc, aes(x=reorder(rownames(sigtab.D4.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 4')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', "#999999"))
deseq.D4.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.fc'
sigtab.D4.fc$OTU <- rownames(sigtab.D4.fc)
sigtab.D4.fc$comp <- 'D4_nasal_IFvsNON'
#If there are duplicate OTU rownames, it'll add an extra "1, 2, 3" etc. in numerical order at the end of the rowname
    
#Create a final table ('final.sigtab') that lists all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D4.fc'. 
#Other within-day comparisons that had significant changes in beta diversity between two treatment groups will be added to 'final.sigtab'.
final.sigtab <- sigtab.D4.fc
    
########## Day 4 Nasal IM vs NON  ####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "4_NW_IM")
#Output:
#IM = 22
sum(meta2$Set == "4_NW_NON")
#Output:
#NON = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.ic' table
FS1b.D4.nw.De
res.D4.ic = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IM","4_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D4.ic = res.D4.ic[which(res.D4.ic$padj < .05), ]
sigtab.D4.ic = cbind(as(sigtab.D4.ic, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.ic), ], "matrix"))
format(sigtab.D4.ic$padj, scientific = TRUE)
sigtab.D4.ic$newp <- format(round(sigtab.D4.ic$padj, digits = 3), scientific = TRUE)
sigtab.D4.ic$Treatment <- ifelse(sigtab.D4.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D4.ic'
deseq.D4.ic <- ggplot(sigtab.D4.ic, aes(x=reorder(rownames(sigtab.D4.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 4')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM","NON"), values = c('#56B4E9', '#999999'))
deseq.D4.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.ic'
sigtab.D4.ic
sigtab.D4.ic$OTU <- rownames(sigtab.D4.ic)
sigtab.D4.ic
sigtab.D4.ic$comp <- 'D4_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D4.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D4.ic)
    
######### Day 4 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "4_NW_IF")
#Output:
#IF = 20
sum(meta2$Set == "4_NW_IM")
#Output:
#IM = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.if' table
res.D4.if = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IM","4_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D4.if
sigtab.D4.if = res.D4.if[which(res.D4.if$padj < .05), ]
sigtab.D4.if = cbind(as(sigtab.D4.if, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.if), ], "matrix"))
format(sigtab.D4.if$padj, scientific = TRUE)
sigtab.D4.if$newp <- format(round(sigtab.D4.if$padj, digits = 3), scientific = TRUE)
sigtab.D4.if$Treatment <- ifelse(sigtab.D4.if$log2FoldChange >=0, "IM", "IF")
sigtab.D4.if
    
#Plot 'sigtab.D4.if'
deseq.D4.if <- ggplot(sigtab.D4.if, aes(x=reorder(rownames(sigtab.D4.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 4')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF' , 'IM'), values = c('#E69F00', '#56B4EF'))
deseq.D4.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.if'
sigtab.D4.if$OTU <- rownames(sigtab.D4.if)
sigtab.D4.if$comp <- 'D4_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D4.if' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D4.if)
    
################################################## Day 7 ############################################################
    
########## Day 7 Nasal #####################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D7.nw.De'
FS1b.D7.nw <- subset_samples(FS1b.genus, Day == 7 & Tissue == 'NW')
sample_sums(FS1b.D7.nw)
FS1b.D7.nw <- prune_taxa(taxa_sums(FS1b.D7.nw) > 1, FS1b.D7.nw)
rowSums(FS1b.D7.nw@otu_table)
FS1b.D7.nw.De <- phyloseq_to_deseq2(FS1b.D7.nw, ~ Set)
FS1b.D7.nw.De <- DESeq(FS1b.D7.nw.De, test = "Wald", fitType = "parametric")
FS1b.D7.nw.De$Set
    
########## Day 7 Nasal IF vs NON #####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "7_NW_IF")
#Output:
#IF = 14
sum(meta2$Set == "7_NW_NON")
#Output:
#NON = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.fc' table
res.D7.fc = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IF","7_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D7.fc = res.D7.fc[which(res.D7.fc$padj < .05), ]
sigtab.D7.fc = cbind(as(sigtab.D7.fc, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.fc), ], "matrix"))
format(sigtab.D7.fc$padj, scientific = TRUE)
sigtab.D7.fc$newp <- format(round(sigtab.D7.fc$padj, digits = 3), scientific = TRUE)
sigtab.D7.fc$Treatment <- ifelse(sigtab.D7.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D7.fc'
deseq.D7.fc <- ggplot(sigtab.D7.fc, aes(x=reorder(rownames(sigtab.D7.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 7')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF","NON"), values = c('#E69F00','#999999'))
deseq.D7.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.fc'
sigtab.D7.fc$OTU <- rownames(sigtab.D7.fc)
sigtab.D7.fc$comp <- 'D7_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D7.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.fc)
    
############# Day 7 Nasal IM vs NON  ######################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "7_NW_IM")
#Output:
#IM = 15
sum(meta2$Set == "7_NW_NON")
#Output:
#NON = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.ic' table
FS1b.D7.nw.De
res.D7.ic = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_NON")) 
sigtab.D7.ic = res.D7.ic[which(res.D7.ic$padj < .05), ]
sigtab.D7.ic = cbind(as(sigtab.D7.ic, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.ic), ], "matrix"))
format(sigtab.D7.ic$padj, scientific = TRUE)
sigtab.D7.ic$newp <- format(round(sigtab.D7.ic$padj, digits = 3), scientific = TRUE)
sigtab.D7.ic$Treatment <- ifelse(sigtab.D7.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D7.ic'
deseq.D7.ic <- ggplot(sigtab.D7.ic, aes(x=reorder(rownames(sigtab.D7.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 7')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM","NON"), values = c('#56B4EF', '#999999'))
deseq.D7.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.ic'
sigtab.D7.ic$OTU <- rownames(sigtab.D7.ic)
sigtab.D7.ic$comp <- 'D7_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D7.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.ic)
    
######### Day 7 Nasal Wash IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "7_NW_IF")
#Output:
#IF = 14
sum(meta2$Set == "7_NW_IM")
#Output:
#IM = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.if' table
res.D7.if = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D7.if
sigtab.D7.if = res.D7.if[which(res.D7.if$padj < .05), ]
sigtab.D7.if = cbind(as(sigtab.D7.if, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.if), ], "matrix"))
format(sigtab.D7.if$padj, scientific = TRUE)
sigtab.D7.if$newp <- format(round(sigtab.D7.if$padj, digits = 3), scientific = TRUE)
sigtab.D7.if$Treatment <- ifelse(sigtab.D7.if$log2FoldChange >=0, "IM", "IF")
sigtab.D7.if
    
#Plot 'sigtab.D7.if'
deseq.D7.if <- ggplot(sigtab.D7.if, aes(x=reorder(rownames(sigtab.D7.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 7')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4EF'))
deseq.D7.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.if'
sigtab.D7.if$OTU <- rownames(sigtab.D7.if)
sigtab.D7.if$comp <- 'D7_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D7.if' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.if)
    
################################################## Day 11 ############################################################
    
############## Day 11 Nasal ###############
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D11.nw.De'
FS1b.D11.nw <- subset_samples(FS1b.genus, Day == 11 & Tissue == 'NW')
sample_sums(FS1b.D11.nw)
FS1b.D11.nw <- prune_taxa(taxa_sums(FS1b.D11.nw) > 1, FS1b.D11.nw)
rowSums(FS1b.D11.nw@otu_table)
FS1b.D11.nw.De <- phyloseq_to_deseq2(FS1b.D11.nw, ~ Set)
FS1b.D11.nw.De <- DESeq(FS1b.D11.nw.De, test = "Wald", fitType = "parametric")
FS1b.D11.nw.De$Set
    
############## Day 11 Nasal IF vs NON ###############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "11_NW_IF")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.fc' table
res.D11.fc = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IF","11_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D11.fc = res.D11.fc[which(res.D11.fc$padj < .05), ]
sigtab.D11.fc = cbind(as(sigtab.D11.fc, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.fc), ], "matrix"))
format(sigtab.D11.fc$padj, scientific = TRUE)
sigtab.D11.fc$newp <- format(round(sigtab.D11.fc$padj, digits = 3), scientific = TRUE)
sigtab.D11.fc$Treatment <- ifelse(sigtab.D11.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D11.fc'
deseq.D11.fc <- ggplot(sigtab.D11.fc, aes(x=reorder(rownames(sigtab.D11.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 11')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
deseq.D11.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.fc'
sigtab.D11.fc$OTU <- rownames(sigtab.D11.fc)
sigtab.D11.fc$comp <- 'D11_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D11.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D11.fc)
    
########### Day 11 Nasal IM vs NON  ############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "11_NW_IM")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.ic' table
FS1b.D11.nw.De$Set
res.D11.ic = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IM","11_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D11.ic = res.D11.ic[which(res.D11.ic$padj < .05), ]
sigtab.D11.ic = cbind(as(sigtab.D11.ic, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.ic), ], "matrix"))
format(sigtab.D11.ic$padj, scientific = TRUE)
sigtab.D11.ic$newp <- format(round(sigtab.D11.ic$padj, digits = 3), scientific = TRUE)
sigtab.D11.ic$Treatment <- ifelse(sigtab.D11.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D11.ic'
deseq.D11.ic <- ggplot(sigtab.D11.ic, aes(x=reorder(rownames(sigtab.D11.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 11')+ coord_flip() +
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4EF', '#999999'))
deseq.D11.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.ic'
sigtab.D11.ic$OTU <- rownames(sigtab.D11.ic)
sigtab.D11.ic$comp <- 'D11_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D11.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D11.ic)
    
######### Day 11 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "11_NW_IM")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_IF")
#Output:
#IM = 7
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.if' table
res.D11.if = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IM","11_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D11.if
sigtab.D11.if = res.D11.if[which(res.D11.if$padj < .05), ]
sigtab.D11.if = cbind(as(sigtab.D11.if, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.if), ], "matrix"))
format(sigtab.D11.if$padj, scientific = TRUE)
sigtab.D11.if$newp <- format(round(sigtab.D11.if$padj, digits = 3), scientific = TRUE)
sigtab.D11.if$Treatment <- ifelse(sigtab.D11.if$log2FoldChange >=0, "IM", "IF")
sigtab.D11.if
    
#Plot 'sigtab.D11.if'
deseq.D11.if <- ggplot(sigtab.D11.if, aes(x=reorder(rownames(sigtab.D11.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 11')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4E9'))
deseq.D11.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.if'
sigtab.D11.if$OTU <- rownames(sigtab.D11.if)
sigtab.D11.if$comp <- 'D11_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D11.if' to 'final.sigtab'
final.sigtab <- rbind(final.sigtab, sigtab.D11.if)
    
################################################## Day 14 ############################################################
    
########### D14 Nasal IF vs NON################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D14.nw.De'
FS1b.D14.nw <- subset_samples(FS1b.genus, Day == 14 & Tissue == 'NW')
sample_sums(FS1b.D14.nw)
FS1b.D14.nw <- prune_taxa(taxa_sums(FS1b.D14.nw) > 1, FS1b.D14.nw)
rowSums(FS1b.D14.nw@otu_table)
FS1b.D14.nw.De <- phyloseq_to_deseq2(FS1b.D14.nw, ~ Set)
FS1b.D14.nw.De <- DESeq(FS1b.D14.nw.De, test = "Wald", fitType = "parametric")
FS1b.D14.nw.De$Set
    
########### D14 Nasal IF vs NON################
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.fc' table
FS1b.D14.nw.De$Set
res.D14.fc = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IF","14_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D14.fc = res.D14.fc[which(res.D14.fc$padj < .05), ]
sigtab.D14.fc = cbind(as(sigtab.D14.fc, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.fc), ], "matrix"))
sigtab.D14.fc
format(sigtab.D14.fc$padj, scientific = TRUE)
sigtab.D14.fc$newp <- format(round(sigtab.D14.fc$padj, digits = 3), scientific = TRUE)
sigtab.D14.fc$Treatment <- ifelse(sigtab.D14.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D14.fc'
deseq.D14.fc <- ggplot(sigtab.D14.fc, aes(x=reorder(rownames(sigtab.D14.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 14')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
deseq.D14.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.fc'
sigtab.D14.fc$OTU <- rownames(sigtab.D14.fc)
sigtab.D14.fc$comp <- 'D14_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D14.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D14.fc)
    
#Write 'final.sigtab' into csv file
write.csv(final.sigtab, file= "Nasal_FinalDiffAbundGenus.csv")
    
########### Day 14 Nasal IM vs NON  ############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "14_NW_IM")
#Output:
#IF = 6
sum(meta2$Set == "14_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.ic' table
FS1b.D14.nw.De$Set
res.D14.ic = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IM","14_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D14.ic = res.D14.ic[which(res.D14.ic$padj < .05), ]
sigtab.D14.ic = cbind(as(sigtab.D14.ic, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.ic), ], "matrix"))
format(sigtab.D14.ic$padj, scientific = TRUE)
sigtab.D14.ic$newp <- format(round(sigtab.D14.ic$padj, digits = 3), scientific = TRUE)
sigtab.D14.ic$Treatment <- ifelse(sigtab.D14.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D14.ic'
deseq.D14.ic <- ggplot(sigtab.D14.ic, aes(x=reorder(rownames(sigtab.D14.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 14')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4E9', '#999999'))
deseq.D14.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.ic'
sigtab.D14.ic$OTU <- rownames(sigtab.D14.ic)
sigtab.D14.ic$comp <- 'D14_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D14.ic' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.ic)
    
######### Day 14 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "14_NW_IM")
#Output:
#IF = 6
sum(meta2$Set == "14_NW_IF")
#Output:
#IM = 7
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.if' table
res.D14.if = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IM","14_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D14.if
sigtab.D14.if = res.D14.if[which(res.D14.if$padj < .05), ]
sigtab.D14.if = cbind(as(sigtab.D14.if, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.if), ], "matrix"))
format(sigtab.D14.if$padj, scientific = TRUE)
sigtab.D14.if$newp <- format(round(sigtab.D14.if$padj, digits = 3), scientific = TRUE)
sigtab.D14.if$Treatment <- ifelse(sigtab.D14.if$log2FoldChange >=0, "IM", "IF")
sigtab.D14.if
    
#Plot 'sigtab.D14.if'
deseq.D14.if <- ggplot(sigtab.D14.if, aes(x=reorder(rownames(sigtab.D14.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF Group at Nasal Site on Day 14')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4E9'))
deseq.D14.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.if'
sigtab.D14.if$OTU <- rownames(sigtab.D14.if)
sigtab.D14.if$comp <- 'D14_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D14.if' to 'final.nonsigtab'
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.if)
    
#Write 'final.nonsigtab' into csv file
write.csv(final.nonsigtab, file= "Nasal_FinalDiffAbundGenus_NonSignificantDays.csv")
    
#######################################################################################################
    
######### Plots of Differentially Abundant Nasal Families and Genera Combined for Each Pairwise Comparison of Treatment Groups ########
    
#IF and NON Log2fold plot Part A
final_fc <- sigtab.D0.fc
final_fc <- rbind(final_fc, sigtab.D4.fc, sigtab.D7.fc, sigtab.D11.fc, sigtab.D14.fc)
final_fc$Family_Genus <- paste(final_fc$Family, final_fc$Genus) #create new column with Family_Genus
fc_plot <- ggplot(final_fc, aes(x=Family_Genus, log2FoldChange, fill = comp)) +
      geom_bar(stat='identity') +
      labs(x="Family Genus", y = "Total log2 Fold Change") +
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=7), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ 
      coord_flip() +
      ggtitle('Differentially Abundant Nasal Wash Families and Genera between IF and NON Groups') + 
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13))
fc_plot
write.csv(final_fc, file= "FS1b_FinalDiffAbundNasalGenus_FC.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_FC.csv in a spreadsheet editor by removing all genera except for Blautia, Lachnospiraceae_unclassified,Roseburia, Lactobacillus, Acinetobacter, Actinobacillus, and Streptococcus.
#Saved modified file as FS1b_FinalDiffAbundNasalGenus_FC_final.csv
    
#IF and NON Log2fold plot Part B
fc <- read.csv('FS1b_FinalDiffAbundNasalGenus_FC_final.csv', header = TRUE, sep = ",")
head(fc[,1:10])
colnames(fc)
fc$DayComp <- sub('_[A-Za-z]+', '\\2', fc$comp)
unique(fc$DayComp)
fc$Day <- sub('_[A-Za-z]+', '\\2', fc$DayComp)
unique(fc$Day) #D4"  "D7"  "D11"
fc$Day = factor(fc$Day, levels=c("D4","D7", "D11"))
levels(fc$Day) #"D4"  "D7"  "D11"
(fc_logfoldplot <- ggplot(data=fc, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, ncol = 2, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 4)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IF"), values = c('#999999', '#E69F00')))
fc_logfoldplot <- fc_logfoldplot + theme(strip.text = element_text(size= 15, face='italic'))
fc_logfoldplot
write.csv(final_fc, file= "FS1b_FinalDiffAbundNasalGenus_FC_final.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_FC_final.csv by removing Roseburia genus and saving as 
#FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv
    
#IF and NON Log2fold Plot Part C
fc1 <- read.csv('FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv', header = TRUE, sep = ",")
head(fc1[,1:10])
colnames(fc1)
fc1$DayComp <- sub('_[A-Za-z]+', '\\2', fc1$comp)
unique(fc1$DayComp)
fc1$Day <- sub('_[A-Za-z]+', '\\2', fc1$DayComp)
unique(fc1$Day) #D4"  "D7"  "D11"
fc1$Day = factor(fc1$Day, levels=c("D4","D7", "D11"))
levels(fc1$Day) #"D4"  "D7"  "D11"
(fc1_logfoldplot <- ggplot(data=fc1, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, ncol = 2, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 3)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IF"), values = c('#999999', '#E69F00')))
fc1_logfoldplot <- fc1_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
fc1_logfoldplot
    
#IM and NON Log2fold Plot Part A
final_ic <- sigtab.D0.ic
final_ic <- rbind(final_ic, sigtab.D4.ic, sigtab.D7.ic, sigtab.D11.ic, sigtab.D14.ic)
final_ic$Family_Genus <- paste(final_ic$Family, final_ic$Genus) #create new column with Family_Genus
ic_plot <- ggplot(final_ic,  aes(x=Family_Genus, log2FoldChange, fill = comp)) +
      geom_bar(stat='identity') +
      labs(x="Family Genus", y = "Total log2 Fold Change") +
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=7), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ 
      coord_flip() +
      ggtitle('Differentially Abundant Nasal Wash Families and Genera between IM and NON Groups') + 
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13))
ic_plot
write.csv(final_ic, file= "FS1b_FinalDiffAbundNasalGenus_IC.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_IC.csv in a spreadsheet editor by removing all genera except for Actinobacillus and Streptococcus.
#Saved modified file as FS1b_FinalDiffAbundNasalGenus_IC_final.csv
    
#IM and NON Log2fold Plot Part B
ic <- read.csv('FS1b_FinalDiffAbundNasalGenus_IC_final.csv', header = TRUE, sep = ",")
head(ic[,1:10])
colnames(ic)
ic$DayComp <- sub('_[A-Za-z]+', '\\2', ic$comp)
unique(ic$DayComp)
ic$Day <- sub('_[A-Za-z]+', '\\2', ic$DayComp)
unique(ic$Day) #"D0"  "D4"  "D7"  "D11" "D14"
ic$Day = factor(ic$Day, levels=c("D0", "D4","D7", "D11", "D14"))
levels(ic$Day) #"D0"  "D4"  "D7"  "D11" "D14"
(ic_logfoldplot <- ggplot(data=ic, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 2)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IM"), values = c('#999999', '#56B4E9')))
ic_logfoldplot <- ic_logfoldplot + theme(strip.text = element_text(size= 15, face='italic'))
ic_logfoldplot
    
#Combine plots 'fc1_logfoldplot' and 'ic_logfoldplot'
ggtwo=plot_grid(fc1_logfoldplot, ic_logfoldplot, labels = c("A", "B"))
ggtwo
    
#Save 'ggtwo' as a .tiff for publication, 500dpi
ggsave("NasalDESeq.tiff", plot=ggtwo, width = 15, height = 5, dpi = 500, units =c("in"))

#################################################################################################################################################################################################

#EIGHTH SECTION: Nasal and Tonsil Microbiota: Genus Abundance
    
#Purpose: This code generates a list of percent total genera found in each treatment group per day for each tissue and creates a bar graph plot of the data 
    
#Files needed:
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#FS1babundsingleton2000metadata.csv

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")
    
otu <- import_mothur(mothur_shared_file = 'FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('FS1babundsingleton2000metadata.csv', header = TRUE, sep = ",")
head(meta)
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

#Calculate the total percent abundance of each genera on each sample (I used JMP to do this) 
#and save results in a spreadsheet editor such as Excel (see D0_Control_NW.genus.xlsx for an example).
#Since we are only interested in genera that are above 2% abundance, 
#calculate total percentage of all other genera that are less than 2% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "Nasal genus.csv" or "Tonsil genus.csv". 
#Create the following columns: Day, Treatment group, Tissue, Percent Abundance, and Genus. 
#Copy the list of genera and their percent abundances from each of the individual Excel files to the respective "Nasal genus.csv" or "Tonsil genus.csv" spreadsheet.
#Fill in the other columns manually (Day, Treatment Group, Tissue). 
#You should have a file similar to "Nasal genus.csv". Continue to the next step.
nasalgen = read.csv("Nasal genus.csv")
tonsilgen = read.csv("Tonsil genus.csv")

#Label nasal genera that have less than 2% abundance as "Other"
nasalgen$More.than.2=as.character(nasalgen$Genus)
str(nasalgen$More.than.2)
nasalgen$More.than.2[nasalgen$Percent.abundance<2]<-"Other"
    
#Rename treatment group names
nasalgen$Treatment2 = nasalgen$Treatment
nasalgen$Treatment2 <- as.character(nasalgen$Treatment2)
nasalgen$Treatment2[nasalgen2$Treatment2 == 'Control'] <- "NON"
nasalgen$Treatment2[nasalgen$Treatment2 == 'Injected'] <- "IM"
nasalgen$Treatment2[nasalgen$Treatment2 == 'Infeed'] <- "IF"
write.csv(nasalgen, file = "Nasal_GenusPercentAbundanceAbove2percent.csv")
    
#To make sure the total percent abundance of all organisms for each day adds up to 100%, 
#modify the percent abundance for "Other" for each day in Nasal_GenusPercentAbundanceAbove2percent.csv in a spreadsheet editor 
#and save as "Nasal_genus.csv"
    
#Create nasal genera plot
nasalgen2 = read.csv("Nasal_genus.csv", header = TRUE)
nasalgen2.plot <- ggplot(data=nasalgen2, aes(x=Treatment2, y=Percent.abundance, fill=Genus)) +
      geom_bar(stat = 'identity') +
      facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
            axis.title.x = element_blank(),
            strip.text = element_text(size= 15),
            axis.text.y = element_text(size=14), 
            axis.title.y = element_text(size=14), 
            legend.text=element_text(size=14), 
            legend.title=element_text(size=14)) +
      scale_fill_igv(name = "Genus") +
      theme(legend.direction = "vertical")
nasalgen2.plot <- nasalgen2.plot + guides(fill= guide_legend(ncol = 1))
nasalgen2.plot <- nasalgen2.plot + theme(legend.text = element_text(face = 'italic'))
nasalgen2.plot
    
#Save 'nw2tt2.combo' as a .tiff for publication, 500dpi
ggsave("NasalTonsilGenera.tiff", plot=nw2tt2.combo, width = 15, height = 7, dpi = 500, units =c("in"))
