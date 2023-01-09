############################################################
#SRD129 16S - Nasal Microbiota: Genus Abundance



#Purpose: Generate a list of percent total genera found in each treatment group per day for each tissue and creates a bar graph plot of the data

#Files needed:
#SRD129Flu.outsingletons.abund.opti_mcc.shared
#SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#SRD129metadata.csv

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)

#######################################################################

#Import files
otu <- import_mothur(mothur_shared_file = './SRD129Flu.outsingletons.abund.opti_mcc.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './SRD129Flu.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
meta <- read.table('./SRD129metadata.csv', header = TRUE, sep = ",")
head(meta)
colnames(meta)[1] <- 'group' #Rename first column of "meta" as "group" temporarily. Will use "group" to set as rownames later and remove the "group" column
meta$Day<- gsub("D", "", meta$Day) #Remove "D"
meta$group <- as.character(meta$group)
head(meta)
phy_meta <- sample_data(meta)
rownames(phy_meta) <- phy_meta$group
head(phy_meta)
phy_meta <- phy_meta[,-1]
head(phy_meta)

#Create phyloseq-class objects with "otu" and "taxo"
IAV <- phyloseq(otu, taxo)
IAV <- merge_phyloseq(IAV, phy_meta)  #This combines the 'phy_meta' metadata with 'IAV' phyloseq object
colnames(tax_table(IAV)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_sums(IAV) #Calculate the sum of all OTUs for each sample,
IAV <- prune_taxa(taxa_sums(IAV) > 2, IAV)  #Removes OTUs that occur less than 2 times globally
IAV.genus <- tax_glom(IAV, 'Genus')
phyla_tab <- as.data.frame(t(IAV.genus@otu_table)) #Transpose 'IAV.genus' by "otu_table"
head(phyla_tab)
IAV.genus@tax_table[,6]
colnames(phyla_tab) <- IAV.genus@tax_table[,6] #Replace column names in phyla_tab from Otuxxxx with Genus names
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #Calculate the proportion of specific phyla per phyla column in 'phyla_tab'
head(phyla_tab2)
phyla_tab2$group <- rownames(phyla_tab2) #Create new column called "group" in 'phyla_tab2' containing rownames
head(phyla_tab2)
fobar <- merge(meta, phyla_tab2, by = 'group') #Merge 'meta' with 'phyla_tab2' by "group"
head(fobar)
fobar.gather <- fobar %>% gather(Genus, value, -(group:Treatment))  #This converts 'fobar' to long-form dataframe.
#This also created new columns "Genus", "value"; it added columns "group" through "Treatment" before "Genus" and "value"
head(fobar.gather)

#Check to see where the extra "group" column is and remove the column
which(colnames(phyla_tab2)=="group") #Results say column 399 is "group" column
phyla_tab3 <- phyla_tab2[,-399] #Drop the 399th column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0] #Keep the columns that have greater than 0 value, QC check

phyla_tab4$group <- rownames(phyla_tab4) #Rename rownames as "group"
fobar2 <- merge(meta, phyla_tab4, by = 'group')
head(fobar2)
fobar2.gather <- fobar2 %>% gather(Genus, value, -(group:Treatment))
head(fobar2.gather)

#Create "All" column with "Day" and "Treatment" in 'fobar2.gather'
fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, sep = '_')

#Count the number of unique items in 'fobar2.gather'. We're interested in the total unique number of genera
fobar2.gather %>% summarise_each(funs(n_distinct)) #113 total unique genera

fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/398))*100) #398 refers to number of Genera
fobar2.gather$percent_tot <- fobar2.gather$value*100
write.csv(fobar2.gather, "./foobar2.gather.csv")



#Subset each "All" group from fobar2.gather and save as individual csv files.
#Calculate the total % percent abundance of each genera on each sample
#calculate total percentage of all other genera that are less than 1% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "nasalgen.csv".

#day 0
D0_Control.genus <- subset(fobar2.gather, All=="0_Control")
write.csv(D0_Control.genus, file="D0_Control.genus.csv")
D0_Control_sum <- data.frame(aggregate(D0_Control.genus$value, by=list(Category=D0_Control.genus$Genus), FUN=sum))
D0_Control_sum$avg <- (as.numeric(D0_Control_sum$x*100/10))
D0_Control_sum$Genus[D0_Control_sum$avg < 1] <- "Other"
D0_Control_sum$Genus[D0_Control_sum$avg >= 1] <- D0_Control_sum$Category[D0_Control_sum$avg >= 1]
D0_Control_sum$Treatment <- c("Control") 
D0_Control_sum$Day <- c("D0")
sum(D0_Control_sum$`x`)
sum(D0_Control_sum$avg)

nasalgen <- D0_Control_sum

D0_IAV.genus <- subset(fobar2.gather, All=="0_IAV")
write.csv(D0_IAV.genus, file="D0_IAV.genus.csv")
D0_IAV_sum <- data.frame(aggregate(D0_IAV.genus$value, by=list(Category=D0_IAV.genus$Genus), FUN=sum))
D0_IAV_sum$avg <- (as.numeric(D0_IAV_sum$x*100/7))
D0_IAV_sum$Genus[D0_IAV_sum$avg < 1] <- "Other"
D0_IAV_sum$Genus[D0_IAV_sum$avg >= 1] <- D0_IAV_sum$Category[D0_IAV_sum$avg >= 1]
D0_IAV_sum$Treatment <- c("IAV") 
D0_IAV_sum$Day <- c("D0")
sum(D0_IAV_sum$`x`)
sum(D0_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D0_IAV_sum)

#day 1
D1_Control.genus <- subset(fobar2.gather, All=="1_Control")
write.csv(D1_Control.genus, file="D1_Control.genus")

D1_Control_sum <- data.frame(aggregate(D1_Control.genus$value, by=list(Category=D1_Control.genus$Genus), FUN=sum))
D1_Control_sum$avg <- (as.numeric(D1_Control_sum$x*100/10))
D1_Control_sum$Genus[D1_Control_sum$avg < 1] <- "Other"
D1_Control_sum$Genus[D1_Control_sum$avg >= 1] <- D1_Control_sum$Category[D1_Control_sum$avg >= 1]
D1_Control_sum$Treatment <- c("Control") 
D1_Control_sum$Day <- c("D1")
sum(D1_Control_sum$`x`)
sum(D1_Control_sum$avg)

nasalgen <- rbind(nasalgen, D1_Control_sum)


D1_IAV.genus <- subset(fobar2.gather, All=="1_IAV")
write.csv(D1_IAV.genus, file="D1_IAV.genus")

D1_IAV_sum <- data.frame(aggregate(D1_IAV.genus$value, by=list(Category=D1_IAV.genus$Genus), FUN=sum))
D1_IAV_sum$avg <- (as.numeric(D1_IAV_sum$x*100/7))
D1_IAV_sum$Genus[D1_IAV_sum$avg < 1] <- "Other"
D1_IAV_sum$Genus[D1_IAV_sum$avg >= 1] <- D1_IAV_sum$Category[D1_IAV_sum$avg >= 1]
D1_IAV_sum$Treatment <- c("IAV") 
D1_IAV_sum$Day <- c("D1")
sum(D1_IAV_sum$`x`)
sum(D1_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D1_IAV_sum)





#day 3
D3_IAV.genus <- subset(fobar2.gather, All=="3_IAV")
write.csv(D3_IAV.genus, file="D3_IAV.genus")

D3_IAV_sum <- data.frame(aggregate(D3_IAV.genus$value, by=list(Category=D3_IAV.genus$Genus), FUN=sum))
D3_IAV_sum$avg <- (as.numeric(D3_IAV_sum$x*100/7))
D3_IAV_sum$Genus[D3_IAV_sum$avg < 1] <- "Other"
D3_IAV_sum$Genus[D3_IAV_sum$avg >= 1] <- D3_IAV_sum$Category[D3_IAV_sum$avg >= 1]
D3_IAV_sum$Treatment <- c("IAV") 
D3_IAV_sum$Day <- c("D3")
sum(D3_IAV_sum$`x`)
sum(D3_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D3_IAV_sum)


D3_Control.genus <- subset(fobar2.gather, All=="3_Control")
write.csv(D3_Control.genus, file="D3_Control.genus")

D3_Control_sum <- data.frame(aggregate(D3_Control.genus$value, by=list(Category=D3_Control.genus$Genus), FUN=sum))
D3_Control_sum$avg <- (as.numeric(D3_Control_sum$x*100/10))
D3_Control_sum$Genus[D3_Control_sum$avg < 1] <- "Other"
D3_Control_sum$Genus[D3_Control_sum$avg >= 1] <- D3_Control_sum$Category[D3_Control_sum$avg >= 1]
D3_Control_sum$Treatment <- c("Control") 
D3_Control_sum$Day <- c("D3")
sum(D3_Control_sum$`x`)
sum(D3_Control_sum$avg)

nasalgen <- rbind(nasalgen, D3_Control_sum)



#day 7
D7_IAV.genus <- subset(fobar2.gather, All=="7_IAV")
write.csv(D7_IAV.genus, file="D7_IAV.genus")

D7_IAV_sum <- data.frame(aggregate(D7_IAV.genus$value, by=list(Category=D7_IAV.genus$Genus), FUN=sum))
D7_IAV_sum$avg <- (as.numeric(D7_IAV_sum$x*100/7))
D7_IAV_sum$Genus[D7_IAV_sum$avg < 1] <- "Other"
D7_IAV_sum$Genus[D7_IAV_sum$avg >= 1] <- D7_IAV_sum$Category[D7_IAV_sum$avg >= 1]
D7_IAV_sum$Treatment <- c("IAV") 
D7_IAV_sum$Day <- c("D7")
sum(D7_IAV_sum$`x`)
sum(D7_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D7_IAV_sum)


D7_Control.genus <- subset(fobar2.gather, All=="7_Control")
write.csv(D7_Control.genus, file="D7_Control.genus")

D7_Control_sum <- data.frame(aggregate(D7_Control.genus$value, by=list(Category=D7_Control.genus$Genus), FUN=sum))
D7_Control_sum$avg <- (as.numeric(D7_Control_sum$x*100/10))
D7_Control_sum$Genus[D7_Control_sum$avg < 1] <- "Other"
D7_Control_sum$Genus[D7_Control_sum$avg >= 1] <- D7_Control_sum$Category[D7_Control_sum$avg >= 1]
D7_Control_sum$Treatment <- c("Control") 
D7_Control_sum$Day <- c("D7")
sum(D7_Control_sum$`x`)
sum(D7_Control_sum$avg)

nasalgen <- rbind(nasalgen, D7_Control_sum)



#day 10
D10_IAV.genus <- subset(fobar2.gather, All=="10_IAV")
write.csv(D10_IAV.genus, file="D10_IAV.genus")

D10_IAV_sum <- data.frame(aggregate(D10_IAV.genus$value, by=list(Category=D10_IAV.genus$Genus), FUN=sum))
D10_IAV_sum$avg <- (as.numeric(D10_IAV_sum$x*100/7))
D10_IAV_sum$Genus[D10_IAV_sum$avg < 1] <- "Other"
D10_IAV_sum$Genus[D10_IAV_sum$avg >= 1] <- D10_IAV_sum$Category[D10_IAV_sum$avg >= 1]
D10_IAV_sum$Treatment <- c("IAV") 
D10_IAV_sum$Day <- c("D10")
sum(D10_IAV_sum$`x`)
sum(D10_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D10_IAV_sum)


D10_Control.genus <- subset(fobar2.gather, All=="10_Control")
write.csv(D10_Control.genus, file="D10_Control.genus")

D10_Control_sum <- data.frame(aggregate(D10_Control.genus$value, by=list(Category=D10_Control.genus$Genus), FUN=sum))
D10_Control_sum$avg <- (as.numeric(D10_Control_sum$x*100/10))
D10_Control_sum$Genus[D10_Control_sum$avg < 1] <- "Other"
D10_Control_sum$Genus[D10_Control_sum$avg >= 1] <- D10_Control_sum$Category[D10_Control_sum$avg >= 1]
D10_Control_sum$Treatment <- c("Control") 
D10_Control_sum$Day <- c("D10")
sum(D10_Control_sum$`x`)
sum(D10_Control_sum$avg)

nasalgen <- rbind(nasalgen, D10_Control_sum)



#day 14
D14_IAV.genus <- subset(fobar2.gather, All=="14_IAV")
write.csv(D14_IAV.genus, file="D14_IAV.genus")

D14_IAV_sum <- data.frame(aggregate(D14_IAV.genus$value, by=list(Category=D14_IAV.genus$Genus), FUN=sum))
D14_IAV_sum$avg <- (as.numeric(D14_IAV_sum$x*100/7))
D14_IAV_sum$Genus[D14_IAV_sum$avg < 1] <- "Other"
D14_IAV_sum$Genus[D14_IAV_sum$avg >= 1] <- D14_IAV_sum$Category[D14_IAV_sum$avg >= 1]
D14_IAV_sum$Treatment <- c("IAV") 
D14_IAV_sum$Day <- c("D14")
sum(D14_IAV_sum$`x`)
sum(D14_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D14_IAV_sum)


D14_Control.genus <- subset(fobar2.gather, All=="14_Control")
write.csv(D14_Control.genus, file="D14_Control.genus")

D14_Control_sum <- data.frame(aggregate(D14_Control.genus$value, by=list(Category=D14_Control.genus$Genus), FUN=sum))
D14_Control_sum$avg <- (as.numeric(D14_Control_sum$x*100/10))
D14_Control_sum$Genus[D14_Control_sum$avg < 1] <- "Other"
D14_Control_sum$Genus[D14_Control_sum$avg >= 1] <- D14_Control_sum$Category[D14_Control_sum$avg >= 1]
D14_Control_sum$Treatment <- c("Control") 
D14_Control_sum$Day <- c("D14")
sum(D14_Control_sum$`x`)
sum(D14_Control_sum$avg)

nasalgen <- rbind(nasalgen, D14_Control_sum)



#day 21
D21_IAV.genus <- subset(fobar2.gather, All=="21_IAV")
write.csv(D21_IAV.genus, file="D21_IAV.genus")

D21_IAV_sum <- data.frame(aggregate(D21_IAV.genus$value, by=list(Category=D21_IAV.genus$Genus), FUN=sum))
D21_IAV_sum$avg <- (as.numeric(D21_IAV_sum$x*100/7))
D21_IAV_sum$Genus[D21_IAV_sum$avg < 1] <- "Other"
D21_IAV_sum$Genus[D21_IAV_sum$avg >= 1] <- D21_IAV_sum$Category[D21_IAV_sum$avg >= 1]
D21_IAV_sum$Treatment <- c("IAV") 
D21_IAV_sum$Day <- c("D21")
sum(D21_IAV_sum$`x`)
sum(D21_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D21_IAV_sum)


D21_Control.genus <- subset(fobar2.gather, All=="21_Control")
write.csv(D21_Control.genus, file="D21_Control.genus")

D21_Control_sum <- data.frame(aggregate(D21_Control.genus$value, by=list(Category=D21_Control.genus$Genus), FUN=sum))
D21_Control_sum$avg <- (as.numeric(D21_Control_sum$x*100/10))
D21_Control_sum$Genus[D21_Control_sum$avg < 1] <- "Other"
D21_Control_sum$Genus[D21_Control_sum$avg >= 1] <- D21_Control_sum$Category[D21_Control_sum$avg >= 1]
D21_Control_sum$Treatment <- c("Control") 
D21_Control_sum$Day <- c("D21")
sum(D21_Control_sum$`x`)
sum(D21_Control_sum$avg)

nasalgen <- rbind(nasalgen, D21_Control_sum)




#day 36
D36_IAV.genus <- subset(fobar2.gather, All=="36_IAV")
write.csv(D36_IAV.genus, file="D36_IAV.genus")

D36_IAV_sum <- data.frame(aggregate(D36_IAV.genus$value, by=list(Category=D36_IAV.genus$Genus), FUN=sum))
D36_IAV_sum$avg <- (as.numeric(D36_IAV_sum$x*100/7))
D36_IAV_sum$Genus[D36_IAV_sum$avg < 1] <- "Other"
D36_IAV_sum$Genus[D36_IAV_sum$avg >= 1] <- D36_IAV_sum$Category[D36_IAV_sum$avg >= 1]
D36_IAV_sum$Treatment <- c("IAV") 
D36_IAV_sum$Day <- c("D36")
sum(D36_IAV_sum$`x`)
sum(D36_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D36_IAV_sum)



D36_Control.genus <- subset(fobar2.gather, All=="36_Control")
write.csv(D36_Control.genus, file="D36_Control.genus")

D36_Control_sum <- data.frame(aggregate(D36_Control.genus$value, by=list(Category=D36_Control.genus$Genus), FUN=sum))
D36_Control_sum$avg <- (as.numeric(D36_Control_sum$x*100/10))
D36_Control_sum$Genus[D36_Control_sum$avg < 1] <- "Other"
D36_Control_sum$Genus[D36_Control_sum$avg >= 1] <- D36_Control_sum$Category[D36_Control_sum$avg >= 1]
D36_Control_sum$Treatment <- c("Control") 
D36_Control_sum$Day <- c("D36")
sum(D36_Control_sum$`x`)
sum(D36_Control_sum$avg)

nasalgen <- rbind(nasalgen, D36_Control_sum)


#day 42
D42_IAV.genus <- subset(fobar2.gather, All=="42_IAV")
write.csv(D42_IAV.genus, file="D42_IAV.genus")

D42_IAV_sum <- data.frame(aggregate(D42_IAV.genus$value, by=list(Category=D42_IAV.genus$Genus), FUN=sum))
D42_IAV_sum$avg <- (as.numeric(D42_IAV_sum$x*100/7))
D42_IAV_sum$Genus[D42_IAV_sum$avg < 1] <- "Other"
D42_IAV_sum$Genus[D42_IAV_sum$avg >= 1] <- D42_IAV_sum$Category[D42_IAV_sum$avg >= 1]
D42_IAV_sum$Treatment <- c("IAV") 
D42_IAV_sum$Day <- c("D42")
sum(D42_IAV_sum$`x`)
sum(D42_IAV_sum$avg)

nasalgen <- rbind(nasalgen, D42_IAV_sum)



D42_Control.genus <- subset(fobar2.gather, All=="42_Control")
write.csv(D42_Control.genus, file="D42_Control.genus")

D42_Control_sum <- data.frame(aggregate(D42_Control.genus$value, by=list(Category=D42_Control.genus$Genus), FUN=sum))
D42_Control_sum$avg <- (as.numeric(D42_Control_sum$x*100/10))
D42_Control_sum$Genus[D42_Control_sum$avg < 1] <- "Other"
D42_Control_sum$Genus[D42_Control_sum$avg >= 1] <- D42_Control_sum$Category[D42_Control_sum$avg >= 1]
D42_Control_sum$Treatment <- c("Control") 
D42_Control_sum$Day <- c("D42")
sum(D42_Control_sum$`x`)
sum(D42_Control_sum$avg)

nasalgen <- rbind(nasalgen, D42_Control_sum)
write.csv(nasalgen, "./nasalgen.csv")



#Create nasal genera plot

## Housekeeping: Order variables and create color schemes
unique(nasalgen$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
nasalgen$Day = factor(nasalgen$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))

GlasbeyScale= c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300", 
                "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C",
                "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", 
                "#FFB79F", "#858567", "#A10300", "#000000")
MokoleScale= c("#696969", "#006400", "#808000", "#483d8b","#b22222", "#008b8b", "#cd853f",
               "#000080", "#32cd32", "#7f007f", "#8fbc8f", "#b03060", "#ff4500",
               "#ffa500", "#ffff00", "#00ff00", "#9400d3", "#00ffff", "#0000ff",
               "#f08080", "#da70d6", "#f0e68c", "#6495ed", "#90ee90", "#ff1493",
               "#7b68ee", "#87cefa")
MokoleScale1= c("#696969", "#006400", "#808000", "#483d8b","#b22222", "#008b8b", "#cd853f",
                "#000080", "#32cd32", "#7f007f", "#8fbc8f", "#b03060", "#ff4500",
                "#ffa500", "#ffff00", "#00ff00", "#9400d3", "#00ffff", "#0000ff",
                "#f08080", "#da70d6", "#f0e68c", "#6495ed", "#90ee90", "#ff1493",
                "#7b68ee", "#152238")


nasalgen$Genus <- factor(nasalgen$Genus, levels = c("Actinobacillus", "Aerococcus", "Blautia", "Carnobacteriaceae_unclassified", 
                                                    "Chryseobacterium", "Clostridium_sensu_stricto_1", "Enhydrobacter", "Lachnospiraceae_unclassified", 
                                                    "Lactobacillus", "Lactococcus", "Moraxella", "Neisseriaceae_unclassified", "Phascolarctobacterium",
                                                    "Prevotella_9", "Prevotellaceae_NK3B31_group", "Prevotellaceae_UCG-003", "Prevotellaceae_unclassified", 
                                                    "Rikenellaceae_RC9_gut_group", "Rothia", "Ruminococcaceae_UCG-005", "Streptococcus", "Succinivibrionaceae_unclassified", 
                                                    "Terrisporobacter", "Treponema_2", "Weeksellaceae_unclassified", "Weissella", "Other"))




###Plot nasalgen
nasalgen$Day = factor(nasalgen$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21",  "D36", "D42"))
nasalgen.3plot <- nasalgen  %>%
  ggplot(aes(x=Treatment, y=avg, fill=Genus)) + geom_col(position = "fill") +
  facet_wrap(~Day, scales = 'free')+
  ggtitle("") +
  ylab('Percent of Total Community') +
  xlab('') + 
  #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
  facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = MokoleScale1, aesthetics = "fill") +
  theme(legend.direction = "vertical") +
  theme(legend.text = element_text(face = 'italic')) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))

ggsave("./FluControl_NasalGeneraScaledTo100OnePercent.tiff", plot=nasalgen.3plot, width = 15, height = 7, dpi = 500, units =c("in"))



# #Nasalgen.2 abundance plot for each day, more than 1% genera. Alternative plot, KMou
# (nasalgen.2plot <- ggplot(data=nasalgen, aes(x=Treatment, y=avg, fill=Genus)) +
#     geom_bar(stat = 'identity') +
#     #geom_bar(stat= 'identity', colour='black') +
#     #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
#     facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(axis.text.x=element_text(angle=45, hjust=1),
#           axis.title.x = element_blank()) +
#     scale_fill_igv(name = "Genus") +
#     theme(legend.direction = "vertical") +
#     theme(legend.text = element_text(face = 'italic')))
# 
# ggsave("./FluControl_NasalGeneraOnePercentNoFilter.tiff", plot=nasalgen.2plot, width = 15, height = 7, dpi = 500, units =c("in"))


# ###Daniel from FSEP1
# fobar2.gather  %>%
#   ggplot(aes(x=Treatment, y=value, fill=Genus)) + geom_col(position = "fill") +
#   facet_wrap(~Day, scales = 'free')+
#   ggtitle("") +
#   ylab('Percent of Total Community') +
#   xlab('') + theme_bw() 
# 
# unique(fobar2.gather$Genus)

