# **SRD129_IAV Project**

This is the repository for scripts and data files pertaining to the research paper "Resilience of swine nasal microbiota to influenza A virus challenge in a longitudinal study".

R version 4.1.2 and accompanying packages were used to run the scripts.<br/>

Fastq files are located in Bioproject PRJNA525911.

### Nomenclature
IAV = influenza A virus (Flu group)

## **Table of contents**
| Chapter | Description |
| -- | -- |
| [data](./data) | Includes data files needed to carry out R analysis |
| [scripts](./scripts) | Text file of commands used for mothur, R scripts for 16S rRNA analysis|

## **Scripts description and the order to run them**
| Order | Script file name | Description |
| -- | -- | -- |
| 1a | 1a_mothur.txt | Process 16S sequence data and generate output for R scripts. |
| 1b | 1b_OTUtable_Flu.R | Create OTU table from mothur output to use for creating phyloseq objects. |
| 2 | 2_phyloseq_Flu.R | Generate phyloseq object to use for 3_alpha_beta_diversity_Flu.R. Run [adonis](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis) function with distance matrices to assess how variation is attributed to different experimental treatments or uncontrolled covariates. |
| 3 | 3_alpha_beta_diversity_Flu.R | Run alpha (Shannon, Inverse Simpson) and beta diversity (generating NMDS, pairwise comparisons) analyses, data visualization. |
| 4 | 4_Nasal_DeSeq2_genus.R | Identify differentially abundant bacterial taxa (genus level) between groups within each day, data visualization. |
| 5 | 5_Nasal_MagnitudeOfChange_Flu.R | Plot the F-statistic from PERMANOVA pairwise comparisons of Control and Flu groups over time. This displays the magnitude of change in the nasal bacterial community structure of the Flu group relative to Control. Also includes scripts for data visualization. |
| 6 | 6_Genus_Flu.R | Generate a list of percent total genera found in each treatment group per day for each tissue, data visualization.  |
