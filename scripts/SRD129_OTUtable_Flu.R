#########################################
#SRD129 OTU table
#Kathy Mou

#NOTES: 
#This code reorganizes the subsampled shared file and taxonomy file from mothur to an OTU table for use with phyloseq

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on desktop or network drive)
#Mac
setwd("~/Desktop/SRD129/SRD129_2000singletons")

#Save image file
save.image(file="SRD129_OTUtable6.29.18.RData")
save.image(file="SRD129_OTUtable_25Jun2020.RData")

#Load image file
load("SRD129_OTUtable6.29.18.RData")
load("SRD129_OTUtable_25Jun2020.RData")
#########################################

######TAXONOMY FILE#####
#Made a copy of SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy file and removed the "Size" column
#Renamed the copy to SRD129Flu.outsingletons.abund.taxonomy.csv
taxonomy <- read.csv("SRD129Flu.outsingletons.abund.taxonomy.csv")
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #remove (100) and variations of that from the 
                                                                 #Taxonomy column
taxonomy[1:6,2] #show Taxonomy column rows 1 through 6
write.csv(taxonomy, file = "SRD129abundsingleton2000taxonomy.csv") #in excel, remove the first column (unnecessary rownumbers)

#####SHARED FILE#####
#Made a copy of SRD129Flu.outsingletons.abund.opti_mcc.0.03.subsample.shared file and removed the "Label" column
#Renamed the copy to SRD129Flu.outsingletons.abund.subsample.shared.csv
shared <- read.csv("SRD129Flu.outsingletons.abund.subsample.shared.csv")
head(shared)
shared[1:6,1]
shared <- t(shared) #transpose
head(shared)
shared[1:6,1]
write.csv(shared, file = 'SRD129abundsingleton2000shared.csv') #in excel, remove the "V*", and "numOtus" rows
#1299 OTUs total

#####SHARED FILE COMBINED WITH TAXONOMY FILE####
shared2 <- read.csv("SRD129abundsingleton2000shared.csv")
colnames(shared2)
colnames(shared2) [1] <- "OTU" #rename first column to OTU
colnames(shared2)
taxonomy2 <- read.csv("SRD129abundsingleton2000taxonomy.csv")
colnames(taxonomy2)
OTUtable <- merge(shared2, taxonomy2, by.x ="OTU", by.y = "OTU") #merge SRD129 shared and taxonomy files via OTU
head(OTUtable)
nrow(OTUtable) #1299
ncol(OTUtable) #200
write.csv(OTUtable, file= "SRD129abundsingleton2000OTUtable.csv") #in excel, remove the first row (numbered rows)
#Copy the columns "OTU", "Mock" and "NTC" columns to a different file
#Delete columns "Mock" and NTC" in SRD129abundsingleton2000OTUtable.csv

#check OTU table
OTUtable2 <-read.csv("SRD129abundsingleton2000OTUtable.csv", stringsAsFactors = FALSE)
head(OTUtable2)
