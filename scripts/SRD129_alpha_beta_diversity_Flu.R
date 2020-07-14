#######################################################################
#SRD129 16S alpha and beta diversity
#Kathy Mou

#Purpose (beta diversity): This code generates non-metric multidimensional scaling ordination based on Bray-Curtis dissimilarities to create NMDS plots, and runs pairwise.adonis function to identify any significant differences in bacterial composition between treatment groups on a given day.

#Purpose: This code calculates alpha diversity metrics (Shannon, Inverse Simpson) that are plotted as box and whisker plots, and uses wilcoxon rank sum test to assess any significant differences in diversity between treatment groups on a given day.

#Files needed:


#NOTES: 
#This code analyzes alpha and beta diversity statistics (and NMDS plots) for nasal samples 
#Flu and Control groups only, no DNEG12, DNEG6
#This script uses files created in "SRD129_phyloseq_each_tissue.R"

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory on Mac desktop:
setwd("~/Desktop/SRD129/SRD129_2000singletons")
setwd("C:/Users/Kathy.Mou/Desktop/SRD129/SRD129_2000singletons")

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
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
library("ggsci")

#Load image file
load("SRD129_alpha_beta_diversity2.RData")
load("SRD129_alpha_beta_diversity_Flu.RData")



#Save image file
save.image(file="SRD129_alpha_beta_diversity2.RData")
save.image(file='SRD129_alpha_beta_diversity_Flu.RData')

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





#veganCovEllipse: no longer in funfuns R package

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



###########################################################################################################


########################################### Setting up 'phyloseqFlu' (nasal no bad samples) into dataframes for NMDS calculation ######################################

flu.sam <- data.frame(phyloseqFlu@sam_data) #make phyloseqFlu sam_data into dataframe
flu.otu <- data.frame(t(phyloseqFlu@otu_table)) #make phyloseqFlu otu_table into dataframe
class(flu.sam) #data.frame
rownames(flu.sam) == row.names(flu.otu) #make sure rownames between flu.sam and flu.otu match exactly. Yes they do.
flu.sam$numOTUS <- rowSums(flu.otu > 1)
head(flu.sam)
dim(flu.sam) #186 6


#Sum of how many samples in each day*treatment group:
nrow(flu.sam[flu.sam$Day == "D0",]) #20
nrow(flu.sam[flu.sam$Day == "D1",]) #20
nrow(flu.sam[flu.sam$Day == "D3",]) #20
nrow(flu.sam[flu.sam$Day == "D7",]) #20
nrow(flu.sam[flu.sam$Day == "D10",]) #20
nrow(flu.sam[flu.sam$Day == "D14",]) #20
nrow(flu.sam[flu.sam$Day == "D21",]) #19
nrow(flu.sam[flu.sam$Day == "D28",]) #9
nrow(flu.sam[flu.sam$Day == "D36",]) #19
nrow(flu.sam[flu.sam$Day == "D42",]) #19



######################---------BETA DIVERSITY STATS AND PLOTS WITH PAIRWISE ADONIS: NASAL---------######################

##NMDS calculation
flu.otu[1:10, 1:10]
flu.NMDS <- NMDS_ellipse(flu.sam, flu.otu, grouping_set = 'All')
#grouping_set: group by "All"
#run the NMDS_ellipse function (at end of script) to add to R functions
#NMDS_ellipse function created by Jules
#also run veganCovEllipse function as well

#Result: [1] "Ordination stress: 0.168876156130503"

#Separating meta data and ellipse data to two lists to make NMDS plot
head(flu.NMDS)
flu.NMDS[[1]]
head(flu.NMDS[[1]])
tail(flu.NMDS[[1]])
flu.NMDS.2 <- flu.NMDS[[1]] 
#flu.NMDS.2 has meta data + MDS calculations. Select this 1st list of flu.NMDS using double brackets
flu.df_ell.2 <- flu.NMDS[[2]]
#flu.df_ell.2 is accessing 2nd list from flu.NMDS that has ellipse calculations
#need two separate lists for making NMDS plot
flu.df_ell.2$group #20 levels
head(flu.df_ell.2)

#Create "Day" and "Treatment" columns within 'flu.df_ell.2' for faceting purposes
flu.df_ell.2$Day <- sub(' [A-Za-z]+', '\\1', flu.df_ell.2$group) #create Day column
#'\\1' Return the first part of the regular expression D[A-z0-9]+ from flu.df_ell.2$group

#Testing regular expressions
#thing1 <- sub('D[0-9]+', '\\1', 'D0_BB')
#thing2 <- sub('_[A-Za-z]+', '\\1', 'D0_BB')
#thing2

flu.df_ell.2$Treatment <- sub('D[A-z0-9]+ ', '\\2', flu.df_ell.2$group) #create Treatment column
#'\\2' Return the second part of the sub expression ([A-Za-z]+) from flu.df_ell.2$group
head(flu.df_ell.2)
tail(flu.df_ell.2)

#Restructure level for flu.NMDS.2 and flu.df_ell.2
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
flu.adon
flu.adon$pairs #list all comparisons in the "pairs" column
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
#use regular expressions
#[A-Za-z] match all capital, lowercase letters
#+ match a whole word and not just one letter (if you didn't have "+")
#c creates the vector, lump all pairs of specific interest groups together
#make vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
goodcomps.flu.adon
goodcomps.flu.adon.2 <- flu.adon[goodcomps.flu.adon,] #rename this vector
goodcomps.flu.adon.2
goodcomps.flu.adon.2$p.adjusted <- p.adjust(goodcomps.flu.adon.2$p.value, method = 'fdr')
goodcomps.flu.adon.2$p.adjusted2 <- round(goodcomps.flu.adon.2$p.adjusted, 3)
goodcomps.flu.adon.2$p.adjusted2[goodcomps.flu.adon.2$p.adjusted2 > 0.05] <- NA
goodcomps.flu.adon.2
write.csv(goodcomps.flu.adon.2, file='goodcomps.flu.adon.csv', row.names=TRUE)


########################---------ALPHA DIVERSITY STATS AND PLOTS WITH PAIRWISE WILCOX: NASAL---------###################

#######################################-----Nasal no bad samples------######################################################

#Calculating alpha diversity metrics: Shannon, Inverse Simpson
flu.sam$shannon <- diversity(otu.meta4) #diversity is vegan function, default index is shannon
#added a shannon index column in flu.sam
flu.sam$invsimpson <- diversity(otu.meta4,index = 'invsimpson')
#use invsimpson since easier to understand than simpson
#don't need to "inverse" simpson values
#simpson: lower the number = higher diversity
flu.sam$numOTUS <- rowSums(otu.meta4 > 1)
levels(sample_data(flu.sam)$Day) #"D0"     "D1"     "D10"    "D14"    "D21"    "D28"    "D3"     "D36"    "D42"    "D7" 
flu.sam$Day = factor(flu.sam$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42")) # Set the level order of values in "Day" column
levels(sample_data(flu.sam)$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"   
head(flu.sam)

#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within nw2.meta3
average.shannon.invsimpson.numOTUs <- aggregate(flu.sam[, 6:8], list(flu.sam$All), mean)
print(average.shannon.invsimpson.numOTUs)

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
write.csv(average.shannon.invsimpson.numOTUs, file="average.shannon.invsimpson.num.OTUs.txt", row.names=TRUE)
flu.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(flu.sam$shannon, flu.sam$All, p.adjust.method = 'none')
print(flu.pairwise.wilcox.shannon.test)
flu.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(flu.sam$invsimpson, flu.sam$All, p.adjust.method = 'none')
print(flu.pairwise.wilcox.invsimpson.test)

#Generate a box and whisker plot of shannon diversity
flu.shannon <- ggplot(data = flu.sam, aes(x=Treatment, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
#free allows each plot to customize scale to the specific data set (no forced scale applied to all plots)
#position = position_dodge2(preserve = 'total') --> fixes the ggplot box width;
#prevents narrow boxes from forming in plot (they become wider)
flu.shannon
#To change treatment group colors
#nw.shan <- nw.shan + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none")


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
#To change treatment group colors
#nw.invsimp.bw <- nw.invsimp.bw + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none")


#combining plots
flu.combalpha <- plot_grid(flu.shannon + theme(legend.position = "none"), flu.invsimp, labels = "AUTO")
flu.combalpha

#Save 'flu.combalpha' as a .tiff for publication, 500dpi
ggsave("Figure_4.tiff", plot=flu.combalpha, width = 10, height = 5, dpi = 500, units =c("in"))
