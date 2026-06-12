setwd("/Users/qcluke/Desktop/Chapter 2/Analysis/Nanopore Sequencing")


### TABLES ---------------------------------------------------------------------------------------
library(tidyverse)
library(vegan)

## illumina 16sv4 
# import otu and tax table
otu_illu_16sv4 = read.delim(file = "Illumina Short-read/illumina_16sv4_otu_table.tsv")
dim(otu_illu_16sv4) 

colnames(otu_illu_16sv4)[1] = "OTU.ID"

tax_illu_16sv4 = read.delim(file = "Illumina Short-read/illumina_16sv4_tax_table.tsv")
dim(tax_illu_16sv4) 

tax_illu_16sv4 = tax_illu_16sv4[-1,]
dim(tax_illu_16sv4) 

colnames(tax_illu_16sv4)[1] = "OTU.ID"

# remove Chloroplast, Mitochondria, and Unassigned
otu_illu_16sv4_cleaned = otu_illu_16sv4[!otu_illu_16sv4$OTU.ID %in% tax_illu_16sv4$OTU.ID[str_detect(tax_illu_16sv4$Taxon,"Chloroplast")], ]
dim(otu_illu_16sv4_cleaned) 
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned[!otu_illu_16sv4_cleaned$OTU.ID %in% tax_illu_16sv4$OTU.ID[str_detect(tax_illu_16sv4$Taxon,"Mitochondria")], ]
dim(otu_illu_16sv4_cleaned) 
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned[!otu_illu_16sv4_cleaned$OTU.ID %in% tax_illu_16sv4$OTU.ID[str_detect(tax_illu_16sv4$Taxon,"Unassigned")], ]
dim(otu_illu_16sv4_cleaned) 

# decontamination
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_C1_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_C2_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_C3_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_T1_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_T2_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_illu_16sv4_cleaned$`F_T3_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_C1_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_C2_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_C3_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_T1_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_T2_NC`, 0)))
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_illu_16sv4_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned %>% select(-contains("_NC"))
dim(otu_illu_16sv4_cleaned) 

# set names
rownames(otu_illu_16sv4_cleaned) = otu_illu_16sv4_cleaned$OTU.ID
otu_illu_16sv4_cleaned = otu_illu_16sv4_cleaned[,-1]
dim(otu_illu_16sv4_cleaned) 

# transform
otu_illu_16sv4_cleaned = t(otu_illu_16sv4_cleaned)
dim(otu_illu_16sv4_cleaned) 

# check reads
rowSums(otu_illu_16sv4_cleaned)

min(rowSums(otu_illu_16sv4_cleaned)) 
mean(rowSums(otu_illu_16sv4_cleaned)) 

# rarefy otu table
otu_illu_16sv4_rarefied = rrarefy(otu_illu_16sv4_cleaned, 40000)
dim(otu_illu_16sv4_rarefied) 
rowSums(otu_illu_16sv4_rarefied)

# remove rare OTUs
otu_illu_16sv4_rarefied = otu_illu_16sv4_rarefied[,which(colSums(otu_illu_16sv4_rarefied)>9)]
rowSums(otu_illu_16sv4_rarefied)
min(colSums(otu_illu_16sv4_rarefied)) 
dim(otu_illu_16sv4_rarefied) 


## illumina 18sv9
otu_illu_18sv9 = read.delim(file = "Illumina Short-read/illumina_18sv9_otu_table.tsv")
dim(otu_illu_18sv9) 

colnames(otu_illu_18sv9)[1] = "OTU.ID"

tax_illu_18sv9 = read.delim(file = "Illumina Short-read/illumina_18sv9_tax_table.tsv")
dim(tax_illu_18sv9) 

tax_illu_18sv9 = tax_illu_18sv9[-1,]
dim(tax_illu_18sv9) 

colnames(tax_illu_18sv9)[1] = "OTU.ID"

# remove Metazoa,Streptophyta,and Unassigned 
otu_illu_18sv9_cleaned = otu_illu_18sv9[!otu_illu_18sv9$OTU.ID %in% tax_illu_18sv9$OTU.ID[str_detect(tax_illu_18sv9$Taxon,"Metazoa")], ]
dim(otu_illu_18sv9_cleaned) 
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned[!otu_illu_18sv9_cleaned$OTU.ID %in% tax_illu_18sv9$OTU.ID[str_detect(tax_illu_18sv9$Taxon,"Streptophyta")], ]
dim(otu_illu_18sv9_cleaned) 
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned[!otu_illu_18sv9_cleaned$OTU.ID %in% tax_illu_18sv9$OTU.ID[str_detect(tax_illu_18sv9$Taxon,"Unassigned")], ]
dim(otu_illu_18sv9_cleaned) 

# decontamination
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_C1_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_C2_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_C3_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_T1_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_T2_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_illu_18sv9_cleaned$`F_T3_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_C1_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_C2_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_C3_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_T1_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_T2_NC`, 0)))
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_illu_18sv9_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned %>% select(-contains("_NC"))
dim(otu_illu_18sv9_cleaned) 

# set names
rownames(otu_illu_18sv9_cleaned) = otu_illu_18sv9_cleaned$OTU.ID
otu_illu_18sv9_cleaned = otu_illu_18sv9_cleaned[,-1]
dim(otu_illu_18sv9_cleaned) 

# transform
otu_illu_18sv9_cleaned = t(otu_illu_18sv9_cleaned)
dim(otu_illu_18sv9_cleaned) 

# check reads
rowSums(otu_illu_18sv9_cleaned)
min(rowSums(otu_illu_18sv9_cleaned)) 
mean(rowSums(otu_illu_18sv9_cleaned)) 

# rarefy otu table
otu_illu_18sv9_rarefied = rrarefy(otu_illu_18sv9_cleaned, 20000)
dim(otu_illu_18sv9_rarefied) 
rowSums(otu_illu_18sv9_rarefied)

# remove rare OTUs
otu_illu_18sv9_rarefied = otu_illu_18sv9_rarefied[,which(colSums(otu_illu_18sv9_rarefied)>9)]
rowSums(otu_illu_18sv9_rarefied)

min(colSums(otu_illu_18sv9_rarefied)) 
dim(otu_illu_18sv9_rarefied) 


## ont 16sv4 
# import otu and tax table
otu_ont_16sv4 = read.delim(file = "ONT Short-read/ont_16sv4_otu_table.tsv")
dim(otu_ont_16sv4)

colnames(otu_ont_16sv4)[1] = "OTU.ID"
colnames(otu_ont_16sv4)[2:49] = colnames(otu_illu_16sv4)[2:49]

tax_ont_16sv4 = read.delim(file = "ONT Short-read/ont_16sv4_tax_table.tsv")
dim(tax_ont_16sv4) 

tax_ont_16sv4 = tax_ont_16sv4[-1,]
dim(tax_ont_16sv4)

colnames(tax_ont_16sv4)[1] = "OTU.ID"

# remove Chloroplast, Mitochondria, and Unassigned
otu_ont_16sv4_cleaned = otu_ont_16sv4[!otu_ont_16sv4$OTU.ID %in% tax_ont_16sv4$OTU.ID[str_detect(tax_ont_16sv4$Taxon,"Chloroplast")], ]
dim(otu_ont_16sv4_cleaned) 
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned[!otu_ont_16sv4_cleaned$OTU.ID %in% tax_ont_16sv4$OTU.ID[str_detect(tax_ont_16sv4$Taxon,"Mitochondria")], ]
dim(otu_ont_16sv4_cleaned) 
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned[!otu_ont_16sv4_cleaned$OTU.ID %in% tax_ont_16sv4$OTU.ID[str_detect(tax_ont_16sv4$Taxon,"Unassigned")], ]
dim(otu_ont_16sv4_cleaned) 

# decontamination
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_C1_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_C2_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_C3_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_T1_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_T2_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_ont_16sv4_cleaned$`F_T3_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_C1_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_C2_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_C3_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_T1_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_T2_NC`, 0)))
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_ont_16sv4_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned %>% select(-contains("_NC"))
dim(otu_ont_16sv4_cleaned) 

# set names
rownames(otu_ont_16sv4_cleaned) = otu_ont_16sv4_cleaned$OTU.ID
otu_ont_16sv4_cleaned = otu_ont_16sv4_cleaned[,-1]
dim(otu_ont_16sv4_cleaned) 

# transform
otu_ont_16sv4_cleaned = t(otu_ont_16sv4_cleaned)
dim(otu_ont_16sv4_cleaned) 

# check reads
rowSums(otu_ont_16sv4_cleaned)
min(rowSums(otu_ont_16sv4_cleaned)) 
mean(rowSums(otu_ont_16sv4_cleaned)) 

# rarefy otu table
otu_ont_16sv4_rarefied = rrarefy(otu_ont_16sv4_cleaned, 40000)
dim(otu_ont_16sv4_rarefied) 
rowSums(otu_ont_16sv4_rarefied)

# remove rare OTUs
otu_ont_16sv4_rarefied = otu_ont_16sv4_rarefied[,which(colSums(otu_ont_16sv4_rarefied)>9)]
rowSums(otu_ont_16sv4_rarefied)
min(colSums(otu_ont_16sv4_rarefied)) 
dim(otu_ont_16sv4_rarefied)


## ont 18sv9
# import otu and tax table
otu_ont_18sv9 = read.delim(file = "ONT Short-read/ont_18sv9_otu_table.tsv")
dim(otu_ont_18sv9) 

colnames(otu_ont_18sv9)[1] = "OTU.ID"
colnames(otu_ont_18sv9)[2:49] = colnames(otu_illu_18sv9)[2:49]

tax_ont_18sv9 = read.delim(file = "ONT Short-read/ont_18sv9_tax_table.tsv")
dim(tax_ont_18sv9) 

tax_ont_18sv9 = tax_ont_18sv9[-1,]
dim(tax_ont_18sv9)

colnames(tax_ont_18sv9)[1] = "OTU.ID"

# remove Metazoa,Streptophyta,and Unassigned 
otu_ont_18sv9_cleaned = otu_ont_18sv9[!otu_ont_18sv9$OTU.ID %in% tax_ont_18sv9$OTU.ID[str_detect(tax_ont_18sv9$Taxon,"Metazoa")], ]
dim(otu_ont_18sv9_cleaned) 
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned[!otu_ont_18sv9_cleaned$OTU.ID %in% tax_ont_18sv9$OTU.ID[str_detect(tax_ont_18sv9$Taxon,"Streptophyta")], ]
dim(otu_ont_18sv9_cleaned) 
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned[!otu_ont_18sv9_cleaned$OTU.ID %in% tax_ont_18sv9$OTU.ID[str_detect(tax_ont_18sv9$Taxon,"Unassigned")], ]
dim(otu_ont_18sv9_cleaned) 

# decontamination
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_C1_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_C2_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_C3_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_T1_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_T2_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_ont_18sv9_cleaned$`F_T3_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_C1_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_C2_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_C3_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_T1_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_T2_NC`, 0)))
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_ont_18sv9_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned %>% select(-contains("_NC"))
dim(otu_ont_18sv9_cleaned) 

# set names
rownames(otu_ont_18sv9_cleaned) = otu_ont_18sv9_cleaned$OTU.ID
otu_ont_18sv9_cleaned = otu_ont_18sv9_cleaned[,-1]
dim(otu_ont_18sv9_cleaned)

# transform
otu_ont_18sv9_cleaned = t(otu_ont_18sv9_cleaned)
dim(otu_ont_18sv9_cleaned)

# check reads
rowSums(otu_ont_18sv9_cleaned)
min(rowSums(otu_ont_18sv9_cleaned)) 
mean(rowSums(otu_ont_18sv9_cleaned)) 

# rarefy otu table
otu_ont_18sv9_rarefied = rrarefy(otu_ont_18sv9_cleaned, 20000)
dim(otu_ont_18sv9_rarefied) 
rowSums(otu_ont_18sv9_rarefied)

# remove rare OTUs
otu_ont_18sv9_rarefied = otu_ont_18sv9_rarefied[,which(colSums(otu_ont_18sv9_rarefied)>9)]
rowSums(otu_ont_18sv9_rarefied)
min(colSums(otu_ont_18sv9_rarefied)) 
dim(otu_ont_18sv9_rarefied) 


### ont fl16s 
otu_ont_fl16s = read.delim(file = "ONT Long-read/ont_fl16s_otu_table.tsv")
dim(otu_ont_fl16s) 

colnames(otu_ont_fl16s)[1] = "OTU.ID"
colnames(otu_ont_fl16s)[2:49] = colnames(otu_illu_16sv4)[2:49]

tax_ont_fl16s = read.delim(file = "ONT Long-read/ont_fl16s_tax_table.tsv")
dim(tax_ont_fl16s) 

tax_ont_fl16s = tax_ont_fl16s[-1,]
dim(tax_ont_fl16s) 

colnames(tax_ont_fl16s)[1] = "OTU.ID"

# remove Chloroplast, Mitochondria, and Unassigned
otu_ont_fl16s_cleaned = otu_ont_fl16s[!otu_ont_fl16s$OTU.ID %in% tax_ont_fl16s$OTU.ID[str_detect(tax_ont_fl16s$Taxon,"Chloroplast")], ]
dim(otu_ont_fl16s_cleaned) 
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned[!otu_ont_fl16s_cleaned$OTU.ID %in% tax_ont_fl16s$OTU.ID[str_detect(tax_ont_fl16s$Taxon,"Mitochondria")], ]
dim(otu_ont_fl16s_cleaned) 
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned[!otu_ont_fl16s_cleaned$OTU.ID %in% tax_ont_fl16s$OTU.ID[str_detect(tax_ont_fl16s$Taxon,"Unassigned")], ]
dim(otu_ont_fl16s_cleaned) 

# decontamination
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_C1_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_C2_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_C3_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_T1_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_T2_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_ont_fl16s_cleaned$`F_T3_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_C1_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_C2_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_C3_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_T1_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_T2_NC`, 0)))
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_ont_fl16s_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned %>% select(-contains("_NC"))
dim(otu_ont_fl16s_cleaned) 

# set names
rownames(otu_ont_fl16s_cleaned) = otu_ont_fl16s_cleaned$OTU.ID
otu_ont_fl16s_cleaned = otu_ont_fl16s_cleaned[,-1]
dim(otu_ont_fl16s_cleaned) 

# transform
otu_ont_fl16s_cleaned = t(otu_ont_fl16s_cleaned)
dim(otu_ont_fl16s_cleaned) #

# check reads
rowSums(otu_ont_fl16s_cleaned)
min(rowSums(otu_ont_fl16s_cleaned)) 
mean(rowSums(otu_ont_fl16s_cleaned)) 

# rarefy otu table
otu_ont_fl16s_rarefied = rrarefy(otu_ont_fl16s_cleaned, 40000)
dim(otu_ont_fl16s_rarefied) 
rowSums(otu_ont_fl16s_rarefied)

# remove rare OTUs
otu_ont_fl16s_rarefied = otu_ont_fl16s_rarefied[,which(colSums(otu_ont_fl16s_rarefied)>9)]
rowSums(otu_ont_fl16s_rarefied)
min(colSums(otu_ont_fl16s_rarefied)) 
dim(otu_ont_fl16s_rarefied) 



## ont fl18s
# import otu and tax table
otu_ont_fl18s = read.delim(file = "ONT Long-read/ont_fl18s_otu_table.tsv")
dim(otu_ont_fl18s) 

colnames(otu_ont_fl18s)[1] = "OTU.ID"
colnames(otu_ont_fl18s)[2:49] = colnames(otu_illu_18sv9)[2:49]

tax_ont_fl18s = read.delim(file = "ONT Long-read/ont_fl18s_tax_table.tsv")
dim(tax_ont_fl18s) 

tax_ont_fl18s = tax_ont_fl18s[-1,]
dim(tax_ont_fl18s) 

colnames(tax_ont_fl18s)[1] = "OTU.ID"

# remove Metazoa,Streptophyta,and Unassigned 
otu_ont_fl18s_cleaned = otu_ont_fl18s[!otu_ont_fl18s$OTU.ID %in% tax_ont_fl18s$OTU.ID[str_detect(tax_ont_fl18s$Taxon,"Metazoa")], ]
dim(otu_ont_fl18s_cleaned) 
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned[!otu_ont_fl18s_cleaned$OTU.ID %in% tax_ont_fl18s$OTU.ID[str_detect(tax_ont_fl18s$Taxon,"Streptophyta")], ]
dim(otu_ont_fl18s_cleaned) 
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned[!otu_ont_fl18s_cleaned$OTU.ID %in% tax_ont_fl18s$OTU.ID[str_detect(tax_ont_fl18s$Taxon,"Unassigned")], ]
dim(otu_ont_fl18s_cleaned) 

# decontamination
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_C1"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_C1_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_C2"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_C2_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_C3"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_C3_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_T1"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_T1_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_T2"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_T2_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("F_T3"), ~ pmax(.-otu_ont_fl18s_cleaned$`F_T3_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_C1"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_C1_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_C2"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_C2_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_C3"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_C3_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_T1"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_T1_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_T2"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_T2_NC`, 0)))
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% 
  mutate(across(starts_with("S_T3"), ~ pmax(.-otu_ont_fl18s_cleaned$`S_T3_NC`, 0)))

# remove NC columns
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned %>% select(-contains("_NC"))
dim(otu_ont_fl18s_cleaned) 

# set names
rownames(otu_ont_fl18s_cleaned) = otu_ont_fl18s_cleaned$OTU.ID
otu_ont_fl18s_cleaned = otu_ont_fl18s_cleaned[,-1]
dim(otu_ont_fl18s_cleaned) 

# transform
otu_ont_fl18s_cleaned = t(otu_ont_fl18s_cleaned)
dim(otu_ont_fl18s_cleaned) 

# check reads
rowSums(otu_ont_fl18s_cleaned)
min(rowSums(otu_ont_fl18s_cleaned))
mean(rowSums(otu_ont_fl18s_cleaned)) 

# rarefy otu table
otu_ont_fl18s_rarefied = rrarefy(otu_ont_fl18s_cleaned, 20000)
dim(otu_ont_fl18s_rarefied) 
rowSums(otu_ont_fl18s_rarefied)

# remove rare OTUs
otu_ont_fl18s_rarefied = otu_ont_fl18s_rarefied[, which(colSums(otu_ont_fl18s_rarefied) > 9)]
rowSums(otu_ont_fl18s_rarefied)
min(colSums(otu_ont_fl18s_rarefied)) 
dim(otu_ont_fl18s_rarefied) 
























### SEQUENCING QUALITY ---------------------------------------------------------------------------------------
library(tidyverse)

## data
quality = read_csv("Quality/Quality.csv")
str(quality)

quality_cleaned = quality %>%
  filter(!str_detect(sample, "NC")) %>%
  mutate(environment = str_sub(sample, 3, 3), 
         site = str_sub(sample, 3, 4), 
         sampling = str_sub(sample, 1, 1), 
         environment_sampling = paste(environment, sampling, sep = "_"),
         method = case_when(type %in% c("16sv4", "18sv9") ~ "short", type %in% c("fl16s", "fl18s") ~ "long", type %in% c("illu16sv4", "illu18sv9") ~ "illu"),
         name = case_when(
           str_starts(sample, "F-C1") ~ paste0("CS1-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "F-C2") ~ paste0("CS2-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "F-C3") ~ paste0("CS3-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "F-T1") ~ paste0("ES1-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "F-T2") ~ paste0("ES2-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "F-T3") ~ paste0("ES3-AF", str_sub(sample, 5, 6)),
           str_starts(sample, "S-C1") ~ paste0("CS1-PS", str_sub(sample, 5, 6)),
           str_starts(sample, "S-C2") ~ paste0("CS2-PS", str_sub(sample, 5, 6)),
           str_starts(sample, "S-C3") ~ paste0("CS3-PS", str_sub(sample, 5, 6)),
           str_starts(sample, "S-T1") ~ paste0("ES1-PS", str_sub(sample, 5, 6)),
           str_starts(sample, "S-T2") ~ paste0("ES2-PS", str_sub(sample, 5, 6)),
           str_starts(sample, "S-T3") ~ paste0("ES3-PS", str_sub(sample, 5, 6))))
quality_cleaned = as.data.frame(quality_cleaned)

quality_cleaned$type = factor(quality_cleaned$type, levels = c("16sv4","18sv9","fl16s","fl18s", "illu16sv4", "illu18sv9"))
quality_cleaned$environment = factor(quality_cleaned$environment, levels = c("C","T"))
quality_cleaned$sampling = factor(quality_cleaned$sampling, levels = c("F","S"))
quality_cleaned$environment_sampling = factor(quality_cleaned$environment_sampling, levels = c("C_F", "C_S", "T_F", "T_S"))
quality_cleaned$method = factor(quality_cleaned$method, levels = c("illu", "short", "long"))
str(quality_cleaned)


## stats
shapiro.test(quality_cleaned$median_quality_score) 
shapiro.test(quality_cleaned$mean_quality_score) 
shapiro.test(quality_cleaned$proportion) 

library(car)
leveneTest(median_quality_score ~ method, data = quality_cleaned)
leveneTest(mean_quality_score ~ method, data = quality_cleaned) 
leveneTest(proportion ~ method, data = quality_cleaned) 

library(dunn.test)
kruskal.test(mean_quality_score ~ method, data = quality_cleaned) 
dunn.test(quality_cleaned$mean_quality_score, quality_cleaned$method, method = "bonferroni")

kruskal.test(proportion ~ method, data = quality_cleaned) 
dunn.test(quality_cleaned$proportion, quality_cleaned$method, method = "bonferroni")

summary(quality_cleaned[quality_cleaned$method == "illu", ]$mean_quality_score)
sd(quality_cleaned[quality_cleaned$method == "illu", ]$mean_quality_score)

summary(quality_cleaned[quality_cleaned$method == "short", ]$mean_quality_score)
sd(quality_cleaned[quality_cleaned$method == "short", ]$mean_quality_score)

summary(quality_cleaned[quality_cleaned$method == "long", ]$mean_quality_score)
sd(quality_cleaned[quality_cleaned$method == "long", ]$mean_quality_score)

summary(quality_cleaned[quality_cleaned$method == "illu", ]$proportion)
sd(quality_cleaned[quality_cleaned$method == "illu", ]$proportion)

summary(quality_cleaned[quality_cleaned$method == "short", ]$proportion)
sd(quality_cleaned[quality_cleaned$method == "short", ]$proportion)

summary(quality_cleaned[quality_cleaned$method == "long", ]$proportion)
sd(quality_cleaned[quality_cleaned$method == "long", ]$proportion)



## plot
library(tidyverse) 
library(ggsignif)

quality_figure_mean = ggplot(quality_cleaned, aes(x = method, y = mean_quality_score)) +
  #geom_hline(yintercept = 20, linetype = "dotted", color = "black", size = 0.3) +
  geom_point(shape = 16, size = 2, position=position_jitter(width = 0.3, height = 0), aes(color = factor(method)), show.legend = F, alpha = 0.5) +
  geom_boxplot(width=0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.5, color = "black") +
  scale_color_manual(values = c("illu" = "#A9A9A9", "short" = "#9BC4E4", "long" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(10, 40), breaks = c(10, 20, 30, 40)) +
  ggtitle("") + xlab("") + ylab("Phred") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x =element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

quality_figure_proportion = ggplot(quality_cleaned, aes(x = method, y = proportion)) +
  #geom_hline(yintercept = 70, linetype = "dotted", color = "black", size = 0.3) +
  geom_point(shape = 16, size = 2, position=position_jitter(width = 0.3, height = 0), aes(color = factor(method)), show.legend = F, alpha = 0.5) +
  geom_boxplot(width=0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.5, color = "black") +
  scale_color_manual(values = c("illu" = "#A9A9A9", "short" = "#9BC4E4", "long" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(60, 100), breaks = c(60, 70, 80, 90, 100)) +
  ggtitle("") + xlab("") + ylab("Proportion") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 





quality_figure_illumina_16sv4_mean = ggplot(quality_cleaned[quality_cleaned$type == "illu16sv4", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

quality_figure_ont_16sv4_mean = ggplot(quality_cleaned[quality_cleaned$type == "16sv4", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("quality_figure_ont_16sv4_mean.pdf", plot = quality_figure_ont_16sv4_mean, width = 5, height = 10, units = "cm", device = "pdf")

quality_figure_ont_fl16s_mean = ggplot(quality_cleaned[quality_cleaned$type == "fl16s", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 





quality_figure_illumina_18sv9_mean = ggplot(quality_cleaned[quality_cleaned$type == "illu18sv9", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

quality_figure_ont_18sv9_mean = ggplot(quality_cleaned[quality_cleaned$type == "18sv9", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

quality_figure_ont_fl18s_mean = ggplot(quality_cleaned[quality_cleaned$type == "fl18s", ], aes(x = mean_quality_score, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = round(mean_quality_score, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Mean Phred score") + 
  ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank())  




quality_figure_illumina_16sv4_proportion = ggplot(quality_cleaned[quality_cleaned$type == "illu16sv4", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = 1) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

quality_figure_ont_16sv4_proportion = ggplot(quality_cleaned[quality_cleaned$type == "16sv4", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

quality_figure_ont_fl16s_proportion = ggplot(quality_cleaned[quality_cleaned$type == "fl16s", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 





quality_figure_illumina_18sv9_proportion = ggplot(quality_cleaned[quality_cleaned$type == "illu18sv9", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = 1) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

quality_figure_ont_18sv9_proportion = ggplot(quality_cleaned[quality_cleaned$type == "18sv9", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

quality_figure_ont_fl18s_proportion = ggplot(quality_cleaned[quality_cleaned$type == "fl18s", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Proportion of *Q20 reads") + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 


## figure
library(patchwork)

quality_figure_separate_mean_16s_combined = quality_figure_illumina_16sv4_mean + quality_figure_ont_16sv4_mean + quality_figure_ont_fl16s_mean + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("quality_figure_separate_mean_16s_combined.pdf", plot = quality_figure_separate_mean_16s_combined, width = 17, height = 10, units = "cm", device = "pdf")

quality_figure_separate_mean_18s_combined = quality_figure_illumina_18sv9_mean + quality_figure_ont_18sv9_mean + quality_figure_ont_fl18s_mean + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("quality_figure_separate_mean_18s_combined.pdf", plot = quality_figure_separate_mean_18s_combined, width = 17, height = 10, units = "cm", device = "pdf")

quality_figure_separate_proportion_16s_combined = quality_figure_illumina_16sv4_proportion + quality_figure_ont_16sv4_proportion + quality_figure_ont_fl16s_proportion + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("quality_figure_separate_proportion_16s_combined.pdf", plot = quality_figure_separate_proportion_16s_combined, width = 17, height = 10, units = "cm", device = "pdf")

quality_figure_separate_proportion_18s_combined = quality_figure_illumina_18sv9_proportion + quality_figure_ont_18sv9_proportion + quality_figure_ont_fl18s_proportion + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("quality_figure_separate_proportion_18s_combined.pdf", plot = quality_figure_separate_proportion_18s_combined, width = 17, height = 10, units = "cm", device = "pdf")

























### TAXONOMIC RESOLUTION -----------------------------------------------------------------------------------------------------------

## data
library(tidyverse)
tax_illu_16sv4_norare = tax_illu_16sv4 %>% filter(OTU.ID %in% colnames(otu_illu_16sv4_rarefied))
dim(otu_illu_16sv4_rarefied) 
dim(tax_illu_16sv4_norare) 

tax_ont_16sv4_norare = tax_ont_16sv4 %>% filter(OTU.ID %in% colnames(otu_ont_16sv4_rarefied))
dim(otu_ont_16sv4_rarefied) 
dim(tax_ont_16sv4_norare) 

tax_ont_fl16s_norare = tax_ont_fl16s %>% filter(OTU.ID %in% colnames(otu_ont_fl16s_rarefied))
dim(otu_ont_fl16s_rarefied) 
dim(tax_ont_fl16s_norare) 

tax_illu_18sv9_norare = tax_illu_18sv9 %>% filter(OTU.ID %in% colnames(otu_illu_18sv9_rarefied))
dim(otu_illu_18sv9_rarefied)
dim(tax_illu_18sv9_norare) 

tax_ont_18sv9_norare = tax_ont_18sv9 %>% filter(OTU.ID %in% colnames(otu_ont_18sv9_rarefied))
dim(otu_ont_18sv9_rarefied) 
dim(tax_ont_18sv9_norare) 

tax_ont_fl18s_norare = tax_ont_fl18s %>% filter(OTU.ID %in% colnames(otu_ont_fl18s_rarefied))
dim(otu_ont_fl18s_rarefied) 
dim(tax_ont_fl18s_norare)

tax = data.frame(
  illu_16sv4_total_otu_number = rowSums(otu_illu_16sv4_rarefied > 0),  
  illu_16sv4_species_level_otu_number = rowSums(otu_illu_16sv4_rarefied[, colnames(otu_illu_16sv4_rarefied) %in% tax_illu_16sv4_norare$OTU.ID[str_detect(tax_illu_16sv4_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0),
  illu_16sv4_species_proportion = rowSums(otu_illu_16sv4_rarefied[, colnames(otu_illu_16sv4_rarefied) %in% tax_illu_16sv4_norare$OTU.ID[str_detect(tax_illu_16sv4_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_illu_16sv4_rarefied > 0) * 100,
  ont_16sv4_total_otu_number = rowSums(otu_ont_16sv4_rarefied > 0),  
  ont_16sv4_species_level_otu_number = rowSums(otu_ont_16sv4_rarefied[, colnames(otu_ont_16sv4_rarefied) %in% tax_ont_16sv4_norare$OTU.ID[str_detect(tax_ont_16sv4_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0),
  ont_16sv4_species_proportion = rowSums(otu_ont_16sv4_rarefied[, colnames(otu_ont_16sv4_rarefied) %in% tax_ont_16sv4_norare$OTU.ID[str_detect(tax_ont_16sv4_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_ont_16sv4_rarefied > 0) * 100, 
  ont_fl16s_total_otu_number = rowSums(otu_ont_fl16s_rarefied > 0),  
  ont_fl16s_species_level_otu_number = rowSums(otu_ont_fl16s_rarefied[, colnames(otu_ont_fl16s_rarefied) %in% tax_ont_fl16s_norare$OTU.ID[str_detect(tax_ont_fl16s_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0),
  ont_fl16s_species_proportion = rowSums(otu_ont_fl16s_rarefied[, colnames(otu_ont_fl16s_rarefied) %in% tax_ont_fl16s_norare$OTU.ID[str_detect(tax_ont_fl16s_norare$Taxon, "s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_ont_fl16s_rarefied > 0) * 100,
  illu_18sv9_total_otu_number = rowSums(otu_illu_18sv9_rarefied > 0),  
  illu_18sv9_species_level_otu_number = rowSums(otu_illu_18sv9_rarefied[, colnames(otu_illu_18sv9_rarefied) %in% tax_illu_18sv9_norare$OTU.ID[str_detect(tax_illu_18sv9_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0),
  illu_18sv9_species_proportion = rowSums(otu_illu_18sv9_rarefied[, colnames(otu_illu_18sv9_rarefied) %in% tax_illu_18sv9_norare$OTU.ID[str_detect(tax_illu_18sv9_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_illu_18sv9_rarefied > 0) * 100,
  ont_18sv9_total_otu_number = rowSums(otu_ont_18sv9_rarefied > 0),  
  ont_18sv9_species_level_otu_number = rowSums(otu_ont_18sv9_rarefied[, colnames(otu_ont_18sv9_rarefied) %in% tax_ont_18sv9_norare$OTU.ID[str_detect(tax_ont_18sv9_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0),
  ont_18sv9_species_proportion = rowSums(otu_ont_18sv9_rarefied[, colnames(otu_ont_18sv9_rarefied) %in% tax_ont_18sv9_norare$OTU.ID[str_detect(tax_ont_18sv9_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_ont_18sv9_rarefied > 0) * 100,
  ont_fl18s_total_otu_number = rowSums(otu_ont_fl18s_rarefied > 0),  
  ont_fl18s_species_level_otu_number = rowSums(otu_ont_fl18s_rarefied[, colnames(otu_ont_fl18s_rarefied) %in% tax_ont_fl18s_norare$OTU.ID[str_detect(tax_ont_fl18s_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0),
  ont_fl18s_species_proportion = rowSums(otu_ont_fl18s_rarefied[, colnames(otu_ont_fl18s_rarefied) %in% tax_ont_fl18s_norare$OTU.ID[str_detect(tax_ont_fl18s_norare$Taxon, "(^|;)s__[^;]*[^\\s;]")]] > 0) / rowSums(otu_ont_fl18s_rarefied > 0) * 100,
)

tax_16s = data.frame(sample = c(rownames(tax), rownames(tax), rownames(tax)), 
                     resolution = c(rep("species", 36), rep("species", 36), rep("species", 36)),
                     proportion = c(tax$illu_16sv4_species_proportion, tax$ont_16sv4_species_proportion, tax$ont_fl16s_species_proportion),
                     method = rep(c("Illumina", "ONT_Short", "ONT_Long"), each = 36))

tax_16s = tax_16s %>% mutate(environment = str_extract(sample, "C|T"), sampling = str_extract(sample, "F|S"), environment_sampling = paste(environment, sampling, sep = "_"))
tax_16s$method = factor (tax_16s$method, levels = c("Illumina","ONT_Short","ONT_Long"))
tax_16s$environment = factor (tax_16s$environment, levels = c("C","T"))
tax_16s$sampling = factor (tax_16s$sampling, levels = c("F","S"))
str(tax_16s) 

tax_16s = tax_16s %>%
  mutate(name = case_when(
           str_starts(sample, "F_C1") ~ paste0("CS1-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "F_C2") ~ paste0("CS2-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "F_C3") ~ paste0("CS3-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "F_T1") ~ paste0("ES1-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "F_T2") ~ paste0("ES2-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "F_T3") ~ paste0("ES3-AF", "-", str_sub(sample, 6)),
           str_starts(sample, "S_C1") ~ paste0("CS1-PS", "-", str_sub(sample, 6)),
           str_starts(sample, "S_C2") ~ paste0("CS2-PS", "-", str_sub(sample, 6)),
           str_starts(sample, "S_C3") ~ paste0("CS3-PS", "-", str_sub(sample, 6)),
           str_starts(sample, "S_T1") ~ paste0("ES1-PS", "-", str_sub(sample, 6)),
           str_starts(sample, "S_T2") ~ paste0("ES2-PS", "-", str_sub(sample, 6)),
           str_starts(sample, "S_T3") ~ paste0("ES3-PS", "-", str_sub(sample, 6))))

str(tax_16s)

tax_18s = data.frame(sample = c(rownames(tax), rownames(tax), rownames(tax)), 
                     resolution = c(rep("species", 36), rep("species", 36), rep("species", 36)),
                     proportion = c(tax$illu_18sv9_species_proportion, tax$ont_18sv9_species_proportion, tax$ont_fl18s_species_proportion),
                     method = rep(c("Illumina", "ONT_Short", "ONT_Long"), each = 36))

tax_18s = tax_18s %>% mutate(environment = str_extract(sample, "C|T"), sampling = str_extract(sample, "F|S"), environment_sampling = paste(environment, sampling, sep = "_"))
tax_18s$method = factor (tax_18s$method, levels = c("Illumina","ONT_Short","ONT_Long"))
tax_18s$environment = factor (tax_18s$environment, levels = c("C","T"))
tax_18s$sampling = factor (tax_18s$sampling, levels = c("F","S"))
str(tax_18s)

tax_18s = tax_18s %>%
  mutate(name = case_when(
    str_starts(sample, "F_C1") ~ paste0("CS1-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "F_C2") ~ paste0("CS2-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "F_C3") ~ paste0("CS3-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "F_T1") ~ paste0("ES1-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "F_T2") ~ paste0("ES2-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "F_T3") ~ paste0("ES3-AF", "-", str_sub(sample, 6)),
    str_starts(sample, "S_C1") ~ paste0("CS1-PS", "-", str_sub(sample, 6)),
    str_starts(sample, "S_C2") ~ paste0("CS2-PS", "-", str_sub(sample, 6)),
    str_starts(sample, "S_C3") ~ paste0("CS3-PS", "-", str_sub(sample, 6)),
    str_starts(sample, "S_T1") ~ paste0("ES1-PS", "-", str_sub(sample, 6)),
    str_starts(sample, "S_T2") ~ paste0("ES2-PS", "-", str_sub(sample, 6)),
    str_starts(sample, "S_T3") ~ paste0("ES3-PS", "-", str_sub(sample, 6))))

str(tax_18s)


## stats
shapiro.test(tax_16s$proportion) 
shapiro.test(tax_18s$proportion) 

library(car)
leveneTest(proportion ~ method, data = tax_16s) 
leveneTest(proportion ~ method, data = tax_18s)

library(dunn.test)
kruskal.test(proportion ~ method, data = tax_16s) 
dunn.test(tax_16s$proportion, tax_16s$method, method = "bonferroni")

kruskal.test(proportion ~ method, data = tax_18s) 
dunn.test(tax_18s$proportion, tax_18s$method, method = "bonferroni")

summary(tax_16s[tax_16s$method == "Illumina", ]$proportion)
sd(tax_16s[tax_16s$method == "Illumina", ]$proportion)
summary(tax_16s[tax_16s$method == "ONT_Short", ]$proportion)
sd(tax_16s[tax_16s$method == "ONT_Short", ]$proportion)
summary(tax_16s[tax_16s$method == "ONT_Long", ]$proportion)
sd(tax_16s[tax_16s$method == "ONT_Long", ]$proportion)

summary(tax_18s[tax_18s$method == "Illumina", ]$proportion)
sd(tax_18s[tax_18s$method == "Illumina", ]$proportion)
summary(tax_18s[tax_18s$method == "ONT_Short", ]$proportion)
sd(tax_18s[tax_18s$method == "ONT_Short", ]$proportion)
summary(tax_18s[tax_18s$method == "ONT_Long", ]$proportion)
sd(tax_18s[tax_18s$method == "ONT_Long", ]$proportion)


##  plot
tax_figure_16s_species_level = ggplot(tax_16s, aes(x = method, y = proportion)) +
  geom_point(size = 2, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.3, height = 0), aes(color = factor(method))) +
  geom_boxplot(width = 0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.6, color = "black") +
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_Short" = "#9BC4E4", "ONT_Long" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(20, 80), breaks = c(20, 40, 60, 80)) +
  ggtitle("") + xlab("") + ylab("Pro") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

tax_figure_18s_species_level = ggplot(tax_18s, aes(x = method, y = proportion)) +
  geom_point(size = 2, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.3, height = 0), aes(color = factor(method))) +
  geom_boxplot(width = 0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.6, color = "black") +
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_Short" = "#9BC4E4", "ONT_Long" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(30, 90), breaks = c(30, 50, 70, 90)) +
  ggtitle("") + xlab("") + ylab("Eu") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))





tax_figure_16s_illu = ggplot(tax_16s[tax_16s$method == "Illumina", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(A) Illumina") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

tax_figure_16s_ontshort = ggplot(tax_16s[tax_16s$method == "ONT_Short", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(B) ONT Short") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

tax_figure_16s_ontlong = ggplot(tax_16s[tax_16s$method == "ONT_Long", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(C) ONT Long") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 



tax_figure_18s_illu = ggplot(tax_18s[tax_18s$method == "Illumina", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(A) Illumina") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

tax_figure_18s_ontshort = ggplot(tax_18s[tax_18s$method == "ONT_Short", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(B) ONT Short") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

tax_figure_18s_ontlong = ggplot(tax_18s[tax_18s$method == "ONT_Long", ], aes(x = proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(proportion, 1)), hjust = -0.2, size = 1.5, family = "Helvetica", fontface = "plain") +
  ggtitle("(C) ONT Long") + 
  xlab("Proportion of species-level OTUs") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 


## figure 
library(patchwork)

figure_tax_res_16s = tax_figure_16s_illu + tax_figure_16s_ontshort + tax_figure_16s_ontlong + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("figure_tax_res_16s.pdf", plot = figure_tax_res_16s, width = 17, height = 10, units = "cm", device = "pdf")

figure_tax_res_18s = tax_figure_18s_illu + tax_figure_18s_ontshort + tax_figure_18s_ontlong + plot_layout(ncol = 3, nrow = 1, byrow = T)

ggsave("figure_tax_res_18s.pdf", plot = figure_tax_res_18s, width = 17, height = 10, units = "cm", device = "pdf")



















### SHARED SPECIES --------------------------------------------------------------------------------
library(tidyverse)
library(ggvenn)

## Data
# define functions
extract_species_name_16s = function(df) {
  df %>%
    mutate(sp_name = str_extract(Taxon, "s__.*")) %>%
    mutate(sp_name = str_remove(sp_name, "s__")) %>%
    mutate(sp_name = str_trim(sp_name)) %>%
    filter(!is.na(sp_name) & sp_name != "") %>%
    pull(sp_name)
}

extract_species_name_18s = function(df) {
  df %>%
    mutate(sp_name = str_extract(Taxon, ";s__.*")) %>%
    mutate(sp_name = str_remove(sp_name, ";s__")) %>%
    mutate(sp_name = str_trim(sp_name)) %>%
    filter(!is.na(sp_name) & sp_name != "") %>%
    pull(sp_name)
}

get_sample_otu_counts_16s <- function(otu_table, tax_table, shared_sp, col_prefix) {
  processed_tax <- tax_table %>%
    mutate(sp_name = str_extract(Taxon, "s__.*")) %>%
    mutate(sp_name = str_remove(sp_name, "s__")) %>%
    mutate(sp_name = str_trim(sp_name)) %>%
    filter(!is.na(sp_name) & sp_name != "")
  all_species_otus <- processed_tax %>%
    filter(!is.na(sp_name) & sp_name != "" & sp_name != "unclassified") %>%
    pull(OTU.ID)
  shared_target_otus <- processed_tax %>%
    filter(sp_name %in% shared_sp) %>%
    pull(OTU.ID)
  valid_all_species <- intersect(all_species_otus, colnames(otu_table))
  valid_shared <- intersect(shared_target_otus, colnames(otu_table))
  presence_matrix <- otu_table > 0
  if(length(valid_all_species) > 0) {
    total_species_present <- rowSums(presence_matrix[, valid_all_species, drop = FALSE])
  } else {
    total_species_present <- rep(0, nrow(otu_table))
  }
  if(length(valid_shared) > 0) {
    shared_otus_present <- rowSums(presence_matrix[, valid_shared, drop = FALSE])
  } else {
    shared_otus_present <- rep(0, nrow(otu_table))
  }
  tibble(
    Sample = rownames(otu_table),
    !!paste0(col_prefix, "_Shared_Count") := shared_otus_present,
    !!paste0(col_prefix, "_Total_Species_Count") := total_species_present, 
    !!paste0(col_prefix, "_Proportion") := shared_otus_present / total_species_present
  )
}

get_sample_otu_counts_18s <- function(otu_table, tax_table, shared_sp, col_prefix) {
  processed_tax <- tax_table %>%
    mutate(sp_name = str_extract(Taxon, ";s__.*")) %>%
    mutate(sp_name = str_remove(sp_name, ";s__")) %>%
    mutate(sp_name = str_trim(sp_name)) %>%
    filter(!is.na(sp_name) & sp_name != "")
  all_species_otus <- processed_tax %>%
    filter(!is.na(sp_name) & sp_name != "" & sp_name != "unclassified") %>%
    pull(OTU.ID)
  shared_target_otus <- processed_tax %>%
    filter(sp_name %in% shared_sp) %>%
    pull(OTU.ID)
  valid_all_species <- intersect(all_species_otus, colnames(otu_table))
  valid_shared <- intersect(shared_target_otus, colnames(otu_table))
  presence_matrix <- otu_table > 0
  if(length(valid_all_species) > 0) {
    total_species_present <- rowSums(presence_matrix[, valid_all_species, drop = FALSE])
  } else {
    total_species_present <- rep(0, nrow(otu_table))
  }
  if(length(valid_shared) > 0) {
    shared_otus_present <- rowSums(presence_matrix[, valid_shared, drop = FALSE])
  } else {
    shared_otus_present <- rep(0, nrow(otu_table))
  }
  tibble(
    Sample = rownames(otu_table),
    !!paste0(col_prefix, "_Shared_Count") := shared_otus_present,
    !!paste0(col_prefix, "_Total_Species_Count") := total_species_present, 
    !!paste0(col_prefix, "_Proportion") := shared_otus_present / total_species_present
  )
}


species_shared_16s_3_methods = intersect(
  intersect(extract_species_name_16s(tax_illu_16sv4_norare), 
            extract_species_name_16s(tax_ont_16sv4_norare)), 
  extract_species_name_16s(tax_ont_fl16s_norare)
)
length(species_shared_16s_3_methods) 


species_shared_16s_proportion_3_methods = get_sample_otu_counts_16s(otu_illu_16sv4_rarefied, tax_illu_16sv4_norare, species_shared_16s_3_methods, "Illu_16SV4") %>%
  full_join(get_sample_otu_counts_16s(otu_ont_16sv4_rarefied, tax_ont_16sv4_norare, species_shared_16s_3_methods, "ONT_16SV4"), by = "Sample") %>%
  full_join(get_sample_otu_counts_16s(otu_ont_fl16s_rarefied, tax_ont_fl16s_norare, species_shared_16s_3_methods, "ONT_FL16S"), by = "Sample")

species_shared_16s_proportion_3_methods = species_shared_16s_proportion_3_methods %>%
  select(Sample, ends_with("_Proportion")) %>%
  pivot_longer(
    cols = -Sample,               
    names_to = "Method",          
    values_to = "Proportion"     
  ) %>%
  mutate(Method = str_remove(Method, "_Proportion")) %>%
  mutate(Method = factor(Method, levels = c("Illu_16SV4", "ONT_16SV4", "ONT_FL16S")),
  Proportion = Proportion * 100)

species_shared_16s_proportion_3_methods = species_shared_16s_proportion_3_methods %>%
  mutate(environment = str_sub(Sample, 3, 3), 
         site = str_sub(Sample, 3, 4), 
         sampling = str_sub(Sample, 1, 1), 
         environment_sampling = paste(environment, sampling, sep = "_"),
         name = case_when(
           str_starts(Sample, "F_C1_") ~ paste0("CS1-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_C2_") ~ paste0("CS2-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_C3_") ~ paste0("CS3-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T1_") ~ paste0("ES1-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T2_") ~ paste0("ES2-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T3_") ~ paste0("ES3-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C1_") ~ paste0("CS1-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C2_") ~ paste0("CS2-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C3_") ~ paste0("CS3-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T1_") ~ paste0("ES1-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T2_") ~ paste0("ES2-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T3_") ~ paste0("ES3-PS-", str_sub(Sample, 6, 6))))
species_shared_16s_proportion_3_methods = as.data.frame(species_shared_16s_proportion_3_methods)
str(species_shared_16s_proportion_3_methods)

species_shared_16s_proportion_3_methods$environment_sampling = factor(species_shared_16s_proportion_3_methods$environment_sampling, levels = c("C_F", "C_S", "T_F", "T_S"))
species_shared_16s_proportion_3_methods$name = factor(species_shared_16s_proportion_3_methods$name, levels = c("CS1-AF-1", "CS1-AF-2", "CS1-AF-3", "CS2-AF-1", "CS2-AF-2", "CS2-AF-3", "CS3-AF-1", "CS3-AF-2", "CS3-AF-3",
                                                               "CS1-PS-1", "CS1-PS-2", "CS1-PS-3", "CS2-PS-1", "CS2-PS-2", "CS2-PS-3", "CS3-PS-1", "CS3-PS-2", "CS3-PS-3",
                                                               "ES1-AF-1", "ES1-AF-2", "ES1-AF-3", "ES2-AF-1", "ES2-AF-2", "ES2-AF-3", "ES3-AF-1", "ES3-AF-2", "ES3-AF-3",
                                                               "ES1-PS-1", "ES1-PS-2", "ES1-PS-3", "ES2-PS-1", "ES2-PS-2", "ES2-PS-3", "ES3-PS-1", "ES3-PS-2", "ES3-PS-3"))
str(species_shared_16s_proportion_3_methods)


species_unique_16s_3_methods = c(setdiff(extract_species_name_16s(tax_illu_16sv4_norare), union(extract_species_name_16s(tax_ont_16sv4_norare), extract_species_name_16s(tax_ont_fl16s_norare))), 
                                 setdiff(extract_species_name_16s(tax_ont_16sv4_norare), union(extract_species_name_16s(tax_illu_16sv4_norare), extract_species_name_16s(tax_ont_fl16s_norare))),
                                 setdiff(extract_species_name_16s(tax_ont_fl16s_norare), union(extract_species_name_16s(tax_illu_16sv4_norare), extract_species_name_16s(tax_ont_16sv4_norare))))




species_shared_18s_3_methods = intersect(
  intersect(extract_species_name_18s(tax_illu_18sv9_norare), 
            extract_species_name_18s(tax_ont_18sv9_norare)), 
  extract_species_name_18s(tax_ont_fl18s_norare)
)
length(species_shared_18s_3_methods) 


species_shared_18s_proportion_3_methods = get_sample_otu_counts_18s(otu_illu_18sv9_rarefied, tax_illu_18sv9_norare, species_shared_18s_3_methods, "Illu_18SV9") %>%
  full_join(get_sample_otu_counts_18s(otu_ont_18sv9_rarefied, tax_ont_18sv9_norare, species_shared_18s_3_methods, "ONT_18SV9"), by = "Sample") %>%
  full_join(get_sample_otu_counts_18s(otu_ont_fl18s_rarefied, tax_ont_fl18s_norare, species_shared_18s_3_methods, "ONT_FL18S"), by = "Sample")

species_shared_18s_proportion_3_methods = species_shared_18s_proportion_3_methods %>%
  select(Sample, ends_with("_Proportion")) %>%
  pivot_longer(
    cols = -Sample,               
    names_to = "Method",          
    values_to = "Proportion"     
  ) %>%
  mutate(Method = str_remove(Method, "_Proportion")) %>%
  mutate(Method = factor(Method, levels = c("Illu_18SV9", "ONT_18SV9", "ONT_FL18S")),
  Proportion = Proportion * 100)

species_shared_18s_proportion_3_methods = species_shared_18s_proportion_3_methods %>%
  mutate(environment = str_sub(Sample, 3, 3), 
         site = str_sub(Sample, 3, 4), 
         sampling = str_sub(Sample, 1, 1), 
         environment_sampling = paste(environment, sampling, sep = "_"),
         name = case_when(
           str_starts(Sample, "F_C1_") ~ paste0("CS1-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_C2_") ~ paste0("CS2-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_C3_") ~ paste0("CS3-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T1_") ~ paste0("ES1-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T2_") ~ paste0("ES2-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "F_T3_") ~ paste0("ES3-AF-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C1_") ~ paste0("CS1-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C2_") ~ paste0("CS2-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_C3_") ~ paste0("CS3-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T1_") ~ paste0("ES1-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T2_") ~ paste0("ES2-PS-", str_sub(Sample, 6, 6)),
           str_starts(Sample, "S_T3_") ~ paste0("ES3-PS-", str_sub(Sample, 6, 6))))
species_shared_18s_proportion_3_methods = as.data.frame(species_shared_18s_proportion_3_methods)
str(species_shared_18s_proportion_3_methods)

species_shared_18s_proportion_3_methods$environment_sampling = factor(species_shared_18s_proportion_3_methods$environment_sampling, levels = c("C_F", "C_S", "T_F", "T_S"))
species_shared_18s_proportion_3_methods$name = factor(species_shared_18s_proportion_3_methods$name, levels = c("CS1-AF-1", "CS1-AF-2", "CS1-AF-3", "CS2-AF-1", "CS2-AF-2", "CS2-AF-3", "CS3-AF-1", "CS3-AF-2", "CS3-AF-3",
                                                                                                               "CS1-PS-1", "CS1-PS-2", "CS1-PS-3", "CS2-PS-1", "CS2-PS-2", "CS2-PS-3", "CS3-PS-1", "CS3-PS-2", "CS3-PS-3",
                                                                                                               "ES1-AF-1", "ES1-AF-2", "ES1-AF-3", "ES2-AF-1", "ES2-AF-2", "ES2-AF-3", "ES3-AF-1", "ES3-AF-2", "ES3-AF-3",
                                                                                                               "ES1-PS-1", "ES1-PS-2", "ES1-PS-3", "ES2-PS-1", "ES2-PS-2", "ES2-PS-3", "ES3-PS-1", "ES3-PS-2", "ES3-PS-3"))
str(species_shared_18s_proportion_3_methods)


species_unique_18s_3_methods = c(setdiff(extract_species_name_18s(tax_illu_18sv9_norare), union(extract_species_name_18s(tax_ont_18sv9_norare), extract_species_name_18s(tax_ont_fl18s_norare))), 
                                 setdiff(extract_species_name_18s(tax_ont_18sv9_norare), union(extract_species_name_18s(tax_illu_18sv9_norare), extract_species_name_18s(tax_ont_fl18s_norare))),
                                 setdiff(extract_species_name_18s(tax_ont_fl18s_norare), union(extract_species_name_18s(tax_illu_18sv9_norare), extract_species_name_18s(tax_ont_18sv9_norare))))


## Stats
shapiro.test(species_shared_16s_proportion_3_methods$Proportion) 
shapiro.test(species_shared_18s_proportion_3_methods$Proportion) 

library(car)
leveneTest(Proportion ~ Method, data = species_shared_16s_proportion_3_methods) 
leveneTest(Proportion ~ Method, data = species_shared_18s_proportion_3_methods) 

library(dunn.test)
kruskal.test(Proportion ~ Method, data = species_shared_16s_proportion_3_methods) 
dunn.test(species_shared_16s_proportion_3_methods$Proportion, species_shared_16s_proportion_3_methods$Method, method = "bonferroni")

kruskal.test(Proportion ~ Method, data = species_shared_18s_proportion_3_methods) 
dunn.test(species_shared_18s_proportion_3_methods$Proportion, species_shared_18s_proportion_3_methods$Method, method = "bonferroni")

summary(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "Illu_16SV4",]$Proportion)
sd(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "Illu_16SV4",]$Proportion)
summary(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_16SV4",]$Proportion)
sd(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_16SV4",]$Proportion)
summary(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_FL16S",]$Proportion)
sd(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_FL16S",]$Proportion)

summary(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "Illu_18SV9",]$Proportion)
sd(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "Illu_18SV9",]$Proportion)
summary(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_18SV9",]$Proportion)
sd(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_18SV9",]$Proportion)
summary(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_FL18S",]$Proportion)
sd(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_FL18S",]$Proportion)


## Plot
species_shared_panel_16s_3_methods = ggplot(species_shared_16s_proportion_3_methods, aes(x = Method, y = Proportion)) +
  geom_point(size = 2, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.3, height = 0), aes(color = factor(Method))) +
  geom_boxplot(width = 0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(Method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.6, color = "black") +
  scale_color_manual(values = c("Illu_16SV4" = "#A9A9A9", "ONT_16SV4" = "#9BC4E4", "ONT_FL16S" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(20, 80), breaks = c(20, 40, 60, 80)) +
  ggtitle("") + xlab("") + ylab("Pro") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))


species_shared_panel_18s_3_methods = ggplot(species_shared_18s_proportion_3_methods, aes(x = Method, y = Proportion)) +
  geom_point(size = 2, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.3, height = 0), aes(color = factor(Method))) +
  geom_boxplot(width = 0.6, fill = NA, outlier.shape = NA, linewidth = 0.6, aes(color = factor(Method)), alpha = 1) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.6, linetype = "solid", linewidth = 0.6, color = "black") +
  scale_color_manual(values = c("Illu_18SV9" = "#A9A9A9", "ONT_18SV9" = "#9BC4E4", "ONT_FL18S" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(30, 90), breaks = c(30, 50, 70, 90)) +
  ggtitle("") + xlab("") + ylab("Eu") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))






species_unique_panel_16s_3_methods = ggplot(species_unique_16s_3_methods_freq, aes(x = reorder(Genus, -Count), y = Count)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", fill = "black", alpha = 1) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  ggtitle("")+ xlab("") + ylab("Pro") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("species_unique_panel_16s_3_methods.pdf", plot = species_unique_panel_16s_3_methods, width = 16, height = 8, units = "cm", device = "pdf")

species_unique_panel_18s_3_methods = ggplot(species_unique_18s_3_methods_freq, aes(x = reorder(Genus, -Count), y = Count)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", fill = "black", alpha = 1) +
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  ggtitle("")+ xlab("") + ylab("Eu") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("species_unique_panel_18s_3_methods.pdf", plot = species_unique_panel_18s_3_methods, width = 16, height = 8, units = "cm", device = "pdf")






species_shared_panel_16s_3_methods_illu = ggplot(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "Illu_16SV4", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

species_shared_panel_16s_3_methods_ontshort = ggplot(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_16SV4", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 

species_shared_panel_16s_3_methods_ontlong = ggplot(species_shared_16s_proportion_3_methods[species_shared_16s_proportion_3_methods$Method == "ONT_FL16S", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 





species_shared_panel_18s_3_methods_illu = ggplot(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "Illu_18SV9", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 4,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

species_shared_panel_18s_3_methods_ontshort = ggplot(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_18SV9", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank())

species_shared_panel_18s_3_methods_ontlong = ggplot(species_shared_18s_proportion_3_methods[species_shared_18s_proportion_3_methods$Method == "ONT_FL18S", ], aes(x = Proportion, y = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, aes(color = factor(environment_sampling), fill = factor(environment_sampling))) +
  scale_color_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_fill_manual(values = c("C_F"= "#8AD9ED", "C_S"= "#05B9E2", "T_F" = "#E0CC9A", "T_S" = "#BB9727")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  geom_text(aes(label = round(Proportion, 1)), size = 1.5, family = "Helvetica", fontface = "plain", hjust = -0.3) +
  ggtitle("") + 
  xlab("Proportion of Common Speicies") + 
  ylab("Sample") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_blank()) 


## Figure 
library(patchwork)

figure_species_shared_16s_3_methods = species_shared_panel_16s_3_methods_illu + species_shared_panel_16s_3_methods_ontshort + species_shared_panel_16s_3_methods_ontlong + plot_layout(ncol = 3, nrow = 1, byrow = T)
ggsave("figure_species_shared_16s_3_methods.pdf", plot = figure_species_shared_16s_3_methods, width = 17, height = 10, units = "cm", device = "pdf")

figure_species_shared_18s_3_methods = species_shared_panel_18s_3_methods_illu + species_shared_panel_18s_3_methods_ontshort + species_shared_panel_18s_3_methods_ontlong + plot_layout(ncol = 3, nrow = 1, byrow = T)
ggsave("figure_species_shared_18s_3_methods.pdf", plot = figure_species_shared_18s_3_methods, width = 17, height = 10, units = "cm", device = "pdf")


figure_quality_resolution_species = quality_figure_mean + quality_figure_proportion +
                                    tax_figure_16s_species_level + tax_figure_18s_species_level +
                                    species_shared_panel_16s_3_methods + species_shared_panel_18s_3_methods + plot_layout(ncol = 2, nrow = 3, byrow = T)
ggsave("figure_quality_resolution_species.pdf", plot = figure_quality_resolution_species, width = 16, height = 18, units = "cm", device = "pdf")



















### DIVERSITY  --------------------------------------------------------------------------------

## Alpha
library(tidyverse)
library(vegan)

# data
alpha = data.frame(OTU_illumina_16sv4 = estimateR(otu_illu_16sv4_rarefied)[1,], 
                   OTU_ont_16sv4 = estimateR(otu_ont_16sv4_rarefied)[1,],
                   OTU_ont_fl16s = estimateR(otu_ont_fl16s_rarefied)[1,],
                   OTU_illumina_18sv9 = estimateR(otu_illu_18sv9_rarefied)[1,], 
                   OTU_ont_18sv9 = estimateR(otu_ont_18sv9_rarefied)[1,],
                   OTU_ont_fl18s = estimateR(otu_ont_fl18s_rarefied)[1,],
                   ACE_illumina_16sv4 = estimateR(otu_illu_16sv4_rarefied)[4,], 
                   ACE_ont_16sv4 = estimateR(otu_ont_16sv4_rarefied)[4,],
                   ACE_ont_fl16s = estimateR(otu_ont_fl16s_rarefied)[4,],
                   ACE_illumina_18sv9 = estimateR(otu_illu_18sv9_rarefied)[4,], 
                   ACE_ont_18sv9 = estimateR(otu_ont_18sv9_rarefied)[4,],
                   ACE_ont_fl18s = estimateR(otu_ont_fl18s_rarefied)[4,],
                   Shannon_illumina_16sv4 = diversity(otu_illu_16sv4_rarefied, index = "shannon"),
                   Shannon_ont_16sv4 = diversity(otu_ont_16sv4_rarefied, index = "shannon"),
                   Shannon_ont_fl16s = diversity(otu_ont_fl16s_rarefied, index = "shannon"),
                   Shannon_illumina_18sv9 = diversity(otu_illu_18sv9_rarefied, index = "shannon"),
                   Shannon_ont_18sv9 = diversity(otu_ont_18sv9_rarefied, index = "shannon"),
                   Shannon_ont_fl18s = diversity(otu_ont_fl18s_rarefied, index = "shannon"),
                   Simpson_illumina_16sv4 = diversity(otu_illu_16sv4_rarefied, index = "simpson"),
                   Simpson_ont_16sv4 = diversity(otu_ont_16sv4_rarefied, index = "simpson"),
                   Simpson_ont_fl16s = diversity(otu_ont_fl16s_rarefied, index = "simpson"),
                   Simpson_illumina_18sv9 = diversity(otu_illu_18sv9_rarefied, index = "simpson"),
                   Simpson_ont_18sv9 = diversity(otu_ont_18sv9_rarefied, index = "simpson"),
                   Simpson_ont_fl18s = diversity(otu_ont_fl18s_rarefied, index = "simpson"),
                   Pielou_illumina_16sv4 = diversity(otu_illu_16sv4_rarefied, index = "shannon")/log(specnumber(otu_illu_16sv4_rarefied)),
                   Pielou_ont_16sv4 = diversity(otu_ont_16sv4_rarefied, index = "shannon")/log(specnumber(otu_ont_16sv4_rarefied)),
                   Pielou_ont_fl16s = diversity(otu_ont_fl16s_rarefied, index = "shannon")/log(specnumber(otu_ont_fl16s_rarefied)),
                   Pielou_illumina_18sv9 = diversity(otu_illu_18sv9_rarefied, index = "shannon")/log(specnumber(otu_illu_18sv9_rarefied)),
                   Pielou_ont_18sv9 = diversity(otu_ont_18sv9_rarefied, index = "shannon")/log(specnumber(otu_ont_18sv9_rarefied)),
                   Pielou_ont_fl18s = diversity(otu_ont_fl18s_rarefied, index = "shannon")/log(specnumber(otu_ont_fl18s_rarefied)))

alpha_OTU_16s = data.frame(sample = c(rownames(alpha), rownames(alpha), rownames(alpha)), 
                           value = c(alpha$OTU_illumina_16sv4, alpha$OTU_ont_16sv4, alpha$OTU_ont_fl16s), 
                           method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 36))
alpha_OTU_16s = alpha_OTU_16s %>% mutate(sampling = str_sub(sample, 1, 1),
                                         environment = str_sub(sample, 3, 3), 
                                         environment_sampling = str_sub(sample, 1, 3))
alpha_OTU_16s$method = factor (alpha_OTU_16s$method, levels = c("Illumina","ONT_short","ONT_long"))
alpha_OTU_16s$sampling = factor (alpha_OTU_16s$sampling, levels = c("F","S"))
alpha_OTU_16s$environment = factor (alpha_OTU_16s$environment, levels = c("C","T"))
str(alpha_OTU_16s)

alpha_OTU_18s = data.frame(sample = c(rownames(alpha), rownames(alpha), rownames(alpha)), 
                           value = c(alpha$OTU_illumina_18sv9, alpha$OTU_ont_18sv9, alpha$OTU_ont_fl18s), 
                           method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 36))
alpha_OTU_18s = alpha_OTU_18s %>% mutate(sampling = str_sub(sample, 1, 1),
                                         environment = str_sub(sample, 3, 3), 
                                         environment_sampling = str_sub(sample, 1, 3))
alpha_OTU_18s$method = factor (alpha_OTU_18s$method, levels = c("Illumina","ONT_short","ONT_long"))
alpha_OTU_18s$sampling = factor (alpha_OTU_18s$sampling, levels = c("F","S"))
alpha_OTU_18s$environment = factor (alpha_OTU_18s$environment, levels = c("C","T"))
str(alpha_OTU_18s)

alpha_shannon_16s = data.frame(sample = c(rownames(alpha), rownames(alpha), rownames(alpha)), 
                               value = c(alpha$Shannon_illumina_16sv4, alpha$Shannon_ont_16sv4, alpha$Shannon_ont_fl16s), 
                               method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 36))
alpha_shannon_16s = alpha_shannon_16s %>% mutate(sampling = str_sub(sample, 1, 1),
                                               environment = str_sub(sample, 3, 3), 
                                               environment_sampling = str_sub(sample, 1, 3))
alpha_shannon_16s$method = factor (alpha_shannon_16s$method, levels = c("Illumina","ONT_short","ONT_long"))
alpha_shannon_16s$sampling = factor (alpha_shannon_16s$sampling, levels = c("F","S"))
alpha_shannon_16s$environment = factor (alpha_shannon_16s$environment, levels = c("C","T"))
str(alpha_shannon_16s)

alpha_shannon_18s = data.frame(sample = c(rownames(alpha), rownames(alpha), rownames(alpha)), 
                               value = c(alpha$Shannon_illumina_18sv9, alpha$Shannon_ont_18sv9, alpha$Shannon_ont_fl18s), 
                               method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 36))
alpha_shannon_18s = alpha_shannon_18s %>% mutate(sampling = str_sub(sample, 1, 1),
                                                 environment = str_sub(sample, 3, 3), 
                                                 environment_sampling = str_sub(sample, 1, 3))
alpha_shannon_18s$method = factor (alpha_shannon_18s$method, levels = c("Illumina","ONT_short","ONT_long"))
alpha_shannon_18s$sampling = factor (alpha_shannon_18s$sampling, levels = c("F","S"))
alpha_shannon_18s$environment = factor (alpha_shannon_18s$environment, levels = c("C","T"))
str(alpha_shannon_18s)


alpha_OTU_16s_statistics = alpha_OTU_16s %>% mutate(location = str_sub(sample, 1, 4)) %>% 
  group_by(method, location, sampling, environment, environment_sampling) %>% 
  summarise(Mean = mean(value), SD = sd(value), SEM = sd(value)/sqrt(3), .groups = "drop") %>% arrange(method, location)

alpha_OTU_16s_statistics$environment = factor (alpha_OTU_16s_statistics$environment, levels = c("T","C"))
str(alpha_OTU_16s_statistics)

alpha_OTU_18s_statistics = alpha_OTU_18s  %>% mutate(location = str_sub(sample, 1, 4)) %>% 
  group_by(method, location, sampling, environment, environment_sampling) %>% 
  summarise(Mean = mean(value), SD = sd(value), SEM = sd(value)/sqrt(3), .groups = "drop") %>% arrange(method, location)

alpha_OTU_18s_statistics$environment = factor (alpha_OTU_18s_statistics$environment, levels = c("T","C"))
str(alpha_OTU_18s_statistics)

alpha_shannon_16s_statistics = alpha_shannon_16s  %>% mutate(location = str_sub(sample, 1, 4)) %>% 
  group_by(method, location, sampling, environment, environment_sampling) %>% 
  summarise(Mean = mean(value), SD = sd(value), SEM = sd(value)/sqrt(3), .groups = "drop") %>% arrange(method, location)

alpha_shannon_16s_statistics$environment = factor (alpha_shannon_16s_statistics$environment, levels = c("T","C"))
str(alpha_shannon_16s_statistics)

alpha_shannon_18s_statistics = alpha_shannon_18s  %>% mutate(location = str_sub(sample, 1, 4)) %>% 
  group_by(method, location, sampling, environment, environment_sampling) %>% 
  summarise(Mean = mean(value), SD = sd(value), SEM = sd(value)/sqrt(3), .groups = "drop") %>% arrange(method, location)

alpha_shannon_18s_statistics$environment = factor (alpha_shannon_18s_statistics$environment, levels = c("T","C"))
str(alpha_shannon_18s_statistics)


# stats
shapiro.test(alpha_OTU_16s_statistics$Mean) 
shapiro.test(alpha_OTU_18s_statistics$Mean)

shapiro.test(alpha_shannon_16s_statistics$Mean) 
shapiro.test(alpha_shannon_18s_statistics$Mean)

library(car)
leveneTest(Mean ~ method, data = alpha_OTU_16s_statistics) 
leveneTest(Mean ~ method, data = alpha_OTU_18s_statistics) 

leveneTest(Mean ~ method, data = alpha_shannon_16s_statistics) 
leveneTest(Mean ~ method, data = alpha_shannon_18s_statistics) 


library(dunn.test)
kruskal.test(Mean ~ method, data = alpha_OTU_16s_statistics) 
dunn.test(alpha_OTU_16s_statistics$Mean, alpha_OTU_16s_statistics$method, method = "bonferroni")

kruskal.test(Mean ~ method, data = alpha_shannon_16s_statistics) 
dunn.test(alpha_shannon_16s_statistics$Mean, alpha_shannon_16s_statistics$method, method = "bonferroni")

kruskal.test(Mean ~ method, data = alpha_OTU_18s_statistics) 
dunn.test(alpha_OTU_18s_statistics$Mean, alpha_OTU_18s_statistics$method, method = "bonferroni")

kruskal.test(Mean ~ method, data = alpha_shannon_18s_statistics) 
dunn.test(alpha_shannon_18s_statistics$Mean, alpha_shannon_18s_statistics$method, method = "bonferroni")


summary(alpha_OTU_16s_statistics[alpha_OTU_16s_statistics$method == "Illumina",]$Mean)
summary(alpha_OTU_16s_statistics[alpha_OTU_16s_statistics$method == "ONT_short",]$Mean)
summary(alpha_OTU_16s_statistics[alpha_OTU_16s_statistics$method == "ONT_long",]$Mean)

summary(alpha_shannon_16s_statistics[alpha_shannon_16s_statistics$method == "Illumina",]$Mean)
summary(alpha_shannon_16s_statistics[alpha_shannon_16s_statistics$method == "ONT_short",]$Mean)
summary(alpha_shannon_16s_statistics[alpha_shannon_16s_statistics$method == "ONT_long",]$Mean)

summary(alpha_OTU_18s_statistics[alpha_OTU_18s_statistics$method == "Illumina",]$Mean)
summary(alpha_OTU_18s_statistics[alpha_OTU_18s_statistics$method == "ONT_short",]$Mean)
summary(alpha_OTU_18s_statistics[alpha_OTU_18s_statistics$method == "ONT_long",]$Mean)

summary(alpha_shannon_18s_statistics[alpha_shannon_18s_statistics$method == "Illumina",]$Mean)
summary(alpha_shannon_18s_statistics[alpha_shannon_18s_statistics$method == "ONT_short",]$Mean)
summary(alpha_shannon_18s_statistics[alpha_shannon_18s_statistics$method == "ONT_long",]$Mean)

# plots
library(tidyverse)
library(ggsignif)

alpha_OTU_16s_average_figure = ggplot(alpha_OTU_16s_statistics, aes(x = method, y = Mean)) +
  geom_point(size = 2, stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 4000), breaks = c(0, 2000, 4000)) +
  ggtitle("") + xlab("") + ylab("Richness") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 


alpha_OTU_18s_average_figure = ggplot(alpha_OTU_18s_statistics, aes(x = method, y = Mean)) +
  geom_point(size = 2, stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 500), breaks = c(0, 250, 500)) +
  ggtitle("") + xlab("") + ylab("Richness") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 


alpha_shannon_16s_average_figure = ggplot(alpha_shannon_16s_statistics, aes(x = method, y = Mean)) +
  geom_point(size = 2, stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) +  
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(1, 7), breaks = c(1, 4, 7)) +
  ggtitle("") + xlab("") + ylab("Shannon") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 


alpha_shannon_18s_average_figure = ggplot(alpha_shannon_18s_statistics, aes(x = method, y = Mean)) +
  geom_point(size = 2, stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(1, 7), breaks = c(1, 4, 7)) +
  ggtitle("") + xlab("") + ylab("Shannon") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 








## beta 
library(tidyverse)
library(vegan)

# data
beta_16s_jaccard_same_sampling = rbind(
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9))
)

row.names(beta_16s_jaccard_same_sampling) = c("ILLUMINA_AF_C1_C2", "ILLUMINA_AF_C1_C3", "ILLUMINA_AF_C1_T1", "ILLUMINA_AF_C1_T2", "ILLUMINA_AF_C1_T3", "ILLUMINA_AF_C2_C3", "ILLUMINA_AF_C2_T1", "ILLUMINA_AF_C2_T2","ILLUMINA_AF_C2_T3", "ILLUMINA_AF_C3_T1", "ILLUMINA_AF_C3_T2", "ILLUMINA_AF_C3_T3", "ILLUMINA_AF_T1_T2", "ILLUMINA_AF_T1_T3","ILLUMINA_AF_T2_T3",
                             "ILLUMINA_PS_C1_C2", "ILLUMINA_PS_C1_C3", "ILLUMINA_PS_C1_T1", "ILLUMINA_PS_C1_T2", "ILLUMINA_PS_C1_T3", "ILLUMINA_PS_C2_C3", "ILLUMINA_PS_C2_T1", "ILLUMINA_PS_C2_T2","ILLUMINA_PS_C2_T3", "ILLUMINA_PS_C3_T1", "ILLUMINA_PS_C3_T2", "ILLUMINA_PS_C3_T3", "ILLUMINA_PS_T1_T2", "ILLUMINA_PS_T1_T3","ILLUMINA_PS_T2_T3",
                             "ONTSHORT_AF_C1_C2", "ONTSHORT_AF_C1_C3", "ONTSHORT_AF_C1_T1", "ONTSHORT_AF_C1_T2", "ONTSHORT_AF_C1_T3", "ONTSHORT_AF_C2_C3", "ONTSHORT_AF_C2_T1", "ONTSHORT_AF_C2_T2","ONTSHORT_AF_C2_T3", "ONTSHORT_AF_C3_T1", "ONTSHORT_AF_C3_T2", "ONTSHORT_AF_C3_T3", "ONTSHORT_AF_T1_T2", "ONTSHORT_AF_T1_T3","ONTSHORT_AF_T2_T3",
                             "ONTSHORT_PS_C1_C2", "ONTSHORT_PS_C1_C3", "ONTSHORT_PS_C1_T1", "ONTSHORT_PS_C1_T2", "ONTSHORT_PS_C1_T3", "ONTSHORT_PS_C2_C3", "ONTSHORT_PS_C2_T1", "ONTSHORT_PS_C2_T2","ONTSHORT_PS_C2_T3", "ONTSHORT_PS_C3_T1", "ONTSHORT_PS_C3_T2", "ONTSHORT_PS_C3_T3", "ONTSHORT_PS_T1_T2", "ONTSHORT_PS_T1_T3","ONTSHORT_PS_T2_T3",
                             "ONTLONG_AF_C1_C2", "ONTLONG_AF_C1_C3", "ONTLONG_AF_C1_T1", "ONTLONG_AF_C1_T2", "ONTLONG_AF_C1_T3", "ONTLONG_AF_C2_C3", "ONTLONG_AF_C2_T1", "ONTLONG_AF_C2_T2","ONTLONG_AF_C2_T3", "ONTLONG_AF_C3_T1", "ONTLONG_AF_C3_T2", "ONTLONG_AF_C3_T3", "ONTLONG_AF_T1_T2", "ONTLONG_AF_T1_T3","ONTLONG_AF_T2_T3",
                             "ONTLONG_PS_C1_C2", "ONTLONG_PS_C1_C3", "ONTLONG_PS_C1_T1", "ONTLONG_PS_C1_T2", "ONTLONG_PS_C1_T3", "ONTLONG_PS_C2_C3", "ONTLONG_PS_C2_T1", "ONTLONG_PS_C2_T2","ONTLONG_PS_C2_T3", "ONTLONG_PS_C3_T1", "ONTLONG_PS_C3_T2", "ONTLONG_PS_C3_T3", "ONTLONG_PS_T1_T2", "ONTLONG_PS_T1_T3","ONTLONG_PS_T2_T3")

beta_16s_jaccard_same_sampling = beta_16s_jaccard_same_sampling %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% mutate(method = str_extract(sample, "ILLUMINA|ONTSHORT|ONTLONG"), sampling = str_extract(sample, "AF|PS"), method_sampling = paste(method, sampling, sep = "_"), between = str_sub(sample, start = -8))
beta_16s_jaccard_same_sampling = beta_16s_jaccard_same_sampling %>% rowwise() %>% mutate(mean = mean(c_across(starts_with("value"))), sd = sd(c_across(starts_with("value"))), sem  = sd / sqrt(9)) %>% ungroup() 

beta_16s_jaccard_same_sampling$method = factor (beta_16s_jaccard_same_sampling$method, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
beta_16s_jaccard_same_sampling$sampling = factor (beta_16s_jaccard_same_sampling$sampling)
beta_16s_jaccard_same_sampling$method_sampling = factor (beta_16s_jaccard_same_sampling$method_sampling)
beta_16s_jaccard_same_sampling$between = factor (beta_16s_jaccard_same_sampling$between)
str(beta_16s_jaccard_same_sampling)

beta_16s_jaccard_same_env = beta_16s_jaccard_same_sampling %>% filter(str_sub(between, 4, 4) == str_sub(between, -2, -2))
beta_16s_jaccard_same_env = beta_16s_jaccard_same_env %>% mutate(environment_sampling = str_sub(between, 2, 4))
beta_16s_jaccard_same_env = beta_16s_jaccard_same_env %>% mutate(environment = str_sub(between, 4, 4))
beta_16s_jaccard_same_env = beta_16s_jaccard_same_env %>%
  pivot_longer(cols = starts_with("value"), 
               names_to = "value_id", 
               values_to = "value") %>%
  select(sample, method, sampling, method_sampling, between, mean, sd, sem, environment_sampling, environment, value_id, value)

beta_16s_jaccard_same_env$environment = factor(beta_16s_jaccard_same_env$environment, levels = c("C", "T"))
str(beta_16s_jaccard_same_env)



beta_16s_bray_same_sampling = rbind(
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_16sv4_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_illu_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_16sv4_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_16sv4_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_16sv4_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl16s_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_fl16s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl16s_rarefied))]), paste0("value", 1:9))
  )

row.names(beta_16s_bray_same_sampling) = c("ILLUMINA_AF_C1_C2", "ILLUMINA_AF_C1_C3", "ILLUMINA_AF_C1_T1", "ILLUMINA_AF_C1_T2", "ILLUMINA_AF_C1_T3", "ILLUMINA_AF_C2_C3", "ILLUMINA_AF_C2_T1", "ILLUMINA_AF_C2_T2","ILLUMINA_AF_C2_T3", "ILLUMINA_AF_C3_T1", "ILLUMINA_AF_C3_T2", "ILLUMINA_AF_C3_T3", "ILLUMINA_AF_T1_T2", "ILLUMINA_AF_T1_T3","ILLUMINA_AF_T2_T3",
                             "ILLUMINA_PS_C1_C2", "ILLUMINA_PS_C1_C3", "ILLUMINA_PS_C1_T1", "ILLUMINA_PS_C1_T2", "ILLUMINA_PS_C1_T3", "ILLUMINA_PS_C2_C3", "ILLUMINA_PS_C2_T1", "ILLUMINA_PS_C2_T2","ILLUMINA_PS_C2_T3", "ILLUMINA_PS_C3_T1", "ILLUMINA_PS_C3_T2", "ILLUMINA_PS_C3_T3", "ILLUMINA_PS_T1_T2", "ILLUMINA_PS_T1_T3","ILLUMINA_PS_T2_T3",
                             "ONTSHORT_AF_C1_C2", "ONTSHORT_AF_C1_C3", "ONTSHORT_AF_C1_T1", "ONTSHORT_AF_C1_T2", "ONTSHORT_AF_C1_T3", "ONTSHORT_AF_C2_C3", "ONTSHORT_AF_C2_T1", "ONTSHORT_AF_C2_T2","ONTSHORT_AF_C2_T3", "ONTSHORT_AF_C3_T1", "ONTSHORT_AF_C3_T2", "ONTSHORT_AF_C3_T3", "ONTSHORT_AF_T1_T2", "ONTSHORT_AF_T1_T3","ONTSHORT_AF_T2_T3",
                             "ONTSHORT_PS_C1_C2", "ONTSHORT_PS_C1_C3", "ONTSHORT_PS_C1_T1", "ONTSHORT_PS_C1_T2", "ONTSHORT_PS_C1_T3", "ONTSHORT_PS_C2_C3", "ONTSHORT_PS_C2_T1", "ONTSHORT_PS_C2_T2","ONTSHORT_PS_C2_T3", "ONTSHORT_PS_C3_T1", "ONTSHORT_PS_C3_T2", "ONTSHORT_PS_C3_T3", "ONTSHORT_PS_T1_T2", "ONTSHORT_PS_T1_T3","ONTSHORT_PS_T2_T3",
                             "ONTLONG_AF_C1_C2", "ONTLONG_AF_C1_C3", "ONTLONG_AF_C1_T1", "ONTLONG_AF_C1_T2", "ONTLONG_AF_C1_T3", "ONTLONG_AF_C2_C3", "ONTLONG_AF_C2_T1", "ONTLONG_AF_C2_T2","ONTLONG_AF_C2_T3", "ONTLONG_AF_C3_T1", "ONTLONG_AF_C3_T2", "ONTLONG_AF_C3_T3", "ONTLONG_AF_T1_T2", "ONTLONG_AF_T1_T3","ONTLONG_AF_T2_T3",
                             "ONTLONG_PS_C1_C2", "ONTLONG_PS_C1_C3", "ONTLONG_PS_C1_T1", "ONTLONG_PS_C1_T2", "ONTLONG_PS_C1_T3", "ONTLONG_PS_C2_C3", "ONTLONG_PS_C2_T1", "ONTLONG_PS_C2_T2","ONTLONG_PS_C2_T3", "ONTLONG_PS_C3_T1", "ONTLONG_PS_C3_T2", "ONTLONG_PS_C3_T3", "ONTLONG_PS_T1_T2", "ONTLONG_PS_T1_T3","ONTLONG_PS_T2_T3")

beta_16s_bray_same_sampling = beta_16s_bray_same_sampling %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% mutate(method = str_extract(sample, "ILLUMINA|ONTSHORT|ONTLONG"), sampling = str_extract(sample, "AF|PS"), method_sampling = paste(method, sampling, sep = "_"), between = str_sub(sample, start = -8))
beta_16s_bray_same_sampling = beta_16s_bray_same_sampling %>% rowwise() %>% mutate(mean = mean(c_across(starts_with("value"))), sd = sd(c_across(starts_with("value"))), sem  = sd / sqrt(9)) %>% ungroup() 

beta_16s_bray_same_sampling$method = factor (beta_16s_bray_same_sampling$method, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
beta_16s_bray_same_sampling$sampling = factor (beta_16s_bray_same_sampling$sampling)
beta_16s_bray_same_sampling$method_sampling = factor (beta_16s_bray_same_sampling$method_sampling)
beta_16s_bray_same_sampling$between = factor (beta_16s_bray_same_sampling$between)
str(beta_16s_bray_same_sampling)

beta_16s_bray_same_env = beta_16s_bray_same_sampling %>% filter(str_sub(between, 4, 4) == str_sub(between, -2, -2))
beta_16s_bray_same_env = beta_16s_bray_same_env %>% mutate(environment_sampling = str_sub(between, 2, 4))
beta_16s_bray_same_env = beta_16s_bray_same_env %>% mutate(environment = str_sub(between, 4, 4))
beta_16s_bray_same_env = beta_16s_bray_same_env %>%
  pivot_longer(cols = starts_with("value"), 
               names_to = "value_id", 
               values_to = "value") %>%
  select(sample, method, sampling, method_sampling, between, mean, sd, sem, environment_sampling, environment, value_id, value)

beta_16s_bray_same_env$environment = factor(beta_16s_bray_same_env$environment, levels = c("C","T"))
str(beta_16s_bray_same_env)


beta_18s_jaccard_same_sampling = rbind(
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "jaccard"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9))
)

row.names(beta_18s_jaccard_same_sampling) = c("ILLUMINA_AF_C1_C2", "ILLUMINA_AF_C1_C3", "ILLUMINA_AF_C1_T1", "ILLUMINA_AF_C1_T2", "ILLUMINA_AF_C1_T3", "ILLUMINA_AF_C2_C3", "ILLUMINA_AF_C2_T1", "ILLUMINA_AF_C2_T2","ILLUMINA_AF_C2_T3", "ILLUMINA_AF_C3_T1", "ILLUMINA_AF_C3_T2", "ILLUMINA_AF_C3_T3", "ILLUMINA_AF_T1_T2", "ILLUMINA_AF_T1_T3","ILLUMINA_AF_T2_T3",
                                "ILLUMINA_PS_C1_C2", "ILLUMINA_PS_C1_C3", "ILLUMINA_PS_C1_T1", "ILLUMINA_PS_C1_T2", "ILLUMINA_PS_C1_T3", "ILLUMINA_PS_C2_C3", "ILLUMINA_PS_C2_T1", "ILLUMINA_PS_C2_T2","ILLUMINA_PS_C2_T3", "ILLUMINA_PS_C3_T1", "ILLUMINA_PS_C3_T2", "ILLUMINA_PS_C3_T3", "ILLUMINA_PS_T1_T2", "ILLUMINA_PS_T1_T3","ILLUMINA_PS_T2_T3",
                                "ONTSHORT_AF_C1_C2", "ONTSHORT_AF_C1_C3", "ONTSHORT_AF_C1_T1", "ONTSHORT_AF_C1_T2", "ONTSHORT_AF_C1_T3", "ONTSHORT_AF_C2_C3", "ONTSHORT_AF_C2_T1", "ONTSHORT_AF_C2_T2","ONTSHORT_AF_C2_T3", "ONTSHORT_AF_C3_T1", "ONTSHORT_AF_C3_T2", "ONTSHORT_AF_C3_T3", "ONTSHORT_AF_T1_T2", "ONTSHORT_AF_T1_T3","ONTSHORT_AF_T2_T3",
                                "ONTSHORT_PS_C1_C2", "ONTSHORT_PS_C1_C3", "ONTSHORT_PS_C1_T1", "ONTSHORT_PS_C1_T2", "ONTSHORT_PS_C1_T3", "ONTSHORT_PS_C2_C3", "ONTSHORT_PS_C2_T1", "ONTSHORT_PS_C2_T2","ONTSHORT_PS_C2_T3", "ONTSHORT_PS_C3_T1", "ONTSHORT_PS_C3_T2", "ONTSHORT_PS_C3_T3", "ONTSHORT_PS_T1_T2", "ONTSHORT_PS_T1_T3","ONTSHORT_PS_T2_T3",
                                "ONTLONG_AF_C1_C2", "ONTLONG_AF_C1_C3", "ONTLONG_AF_C1_T1", "ONTLONG_AF_C1_T2", "ONTLONG_AF_C1_T3", "ONTLONG_AF_C2_C3", "ONTLONG_AF_C2_T1", "ONTLONG_AF_C2_T2","ONTLONG_AF_C2_T3", "ONTLONG_AF_C3_T1", "ONTLONG_AF_C3_T2", "ONTLONG_AF_C3_T3", "ONTLONG_AF_T1_T2", "ONTLONG_AF_T1_T3","ONTLONG_AF_T2_T3",
                                "ONTLONG_PS_C1_C2", "ONTLONG_PS_C1_C3", "ONTLONG_PS_C1_T1", "ONTLONG_PS_C1_T2", "ONTLONG_PS_C1_T3", "ONTLONG_PS_C2_C3", "ONTLONG_PS_C2_T1", "ONTLONG_PS_C2_T2","ONTLONG_PS_C2_T3", "ONTLONG_PS_C3_T1", "ONTLONG_PS_C3_T2", "ONTLONG_PS_C3_T3", "ONTLONG_PS_T1_T2", "ONTLONG_PS_T1_T3","ONTLONG_PS_T2_T3")

beta_18s_jaccard_same_sampling = beta_18s_jaccard_same_sampling %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% mutate(method = str_extract(sample, "ILLUMINA|ONTSHORT|ONTLONG"), sampling = str_extract(sample, "AF|PS"), method_sampling = paste(method, sampling, sep = "_"), between = str_sub(sample, start = -8))
beta_18s_jaccard_same_sampling = beta_18s_jaccard_same_sampling %>% rowwise() %>% mutate(mean = mean(c_across(starts_with("value"))), sd = sd(c_across(starts_with("value"))), sem  = sd / sqrt(9)) %>% ungroup() 

beta_18s_jaccard_same_sampling$method = factor (beta_18s_jaccard_same_sampling$method, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
beta_18s_jaccard_same_sampling$sampling = factor (beta_18s_jaccard_same_sampling$sampling)
beta_18s_jaccard_same_sampling$method_sampling = factor (beta_18s_jaccard_same_sampling$method_sampling)
beta_18s_jaccard_same_sampling$between = factor (beta_18s_jaccard_same_sampling$between)
str(beta_18s_jaccard_same_sampling)

beta_18s_jaccard_same_env = beta_18s_jaccard_same_sampling %>% filter(str_sub(between, 4, 4) == str_sub(between, -2, -2))
beta_18s_jaccard_same_env = beta_18s_jaccard_same_env %>% mutate(environment_sampling = str_sub(between, 2, 4))
beta_18s_jaccard_same_env = beta_18s_jaccard_same_env %>% mutate(environment = str_sub(between, 4, 4))
beta_18s_jaccard_same_env = beta_18s_jaccard_same_env %>%
  pivot_longer(cols = starts_with("value"), 
               names_to = "value_id", 
               values_to = "value") %>%
  select(sample, method, sampling, method_sampling, between, mean, sd, sem, environment_sampling, environment, value_id, value)

beta_18s_jaccard_same_env$environment = factor(beta_18s_jaccard_same_env$environment, levels = c("C","T"))
str(beta_18s_jaccard_same_env)



beta_18s_bray_same_sampling = rbind(
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_illu_18sv9_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_illu_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_illu_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_18sv9_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_18sv9_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_18sv9_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^F_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^F_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_C3_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_T1_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9)),
  setNames(c(as.matrix(vegdist(otu_ont_fl18s_rarefied, "bray"))[grepl("^S_T2_[1-3]$", rownames(otu_ont_fl18s_rarefied)), grepl("^S_T3_[1-3]$", rownames(otu_ont_fl18s_rarefied))]), paste0("value", 1:9))
)

row.names(beta_18s_bray_same_sampling) = c("ILLUMINA_AF_C1_C2", "ILLUMINA_AF_C1_C3", "ILLUMINA_AF_C1_T1", "ILLUMINA_AF_C1_T2", "ILLUMINA_AF_C1_T3", "ILLUMINA_AF_C2_C3", "ILLUMINA_AF_C2_T1", "ILLUMINA_AF_C2_T2","ILLUMINA_AF_C2_T3", "ILLUMINA_AF_C3_T1", "ILLUMINA_AF_C3_T2", "ILLUMINA_AF_C3_T3", "ILLUMINA_AF_T1_T2", "ILLUMINA_AF_T1_T3","ILLUMINA_AF_T2_T3",
                             "ILLUMINA_PS_C1_C2", "ILLUMINA_PS_C1_C3", "ILLUMINA_PS_C1_T1", "ILLUMINA_PS_C1_T2", "ILLUMINA_PS_C1_T3", "ILLUMINA_PS_C2_C3", "ILLUMINA_PS_C2_T1", "ILLUMINA_PS_C2_T2","ILLUMINA_PS_C2_T3", "ILLUMINA_PS_C3_T1", "ILLUMINA_PS_C3_T2", "ILLUMINA_PS_C3_T3", "ILLUMINA_PS_T1_T2", "ILLUMINA_PS_T1_T3","ILLUMINA_PS_T2_T3",
                             "ONTSHORT_AF_C1_C2", "ONTSHORT_AF_C1_C3", "ONTSHORT_AF_C1_T1", "ONTSHORT_AF_C1_T2", "ONTSHORT_AF_C1_T3", "ONTSHORT_AF_C2_C3", "ONTSHORT_AF_C2_T1", "ONTSHORT_AF_C2_T2","ONTSHORT_AF_C2_T3", "ONTSHORT_AF_C3_T1", "ONTSHORT_AF_C3_T2", "ONTSHORT_AF_C3_T3", "ONTSHORT_AF_T1_T2", "ONTSHORT_AF_T1_T3","ONTSHORT_AF_T2_T3",
                             "ONTSHORT_PS_C1_C2", "ONTSHORT_PS_C1_C3", "ONTSHORT_PS_C1_T1", "ONTSHORT_PS_C1_T2", "ONTSHORT_PS_C1_T3", "ONTSHORT_PS_C2_C3", "ONTSHORT_PS_C2_T1", "ONTSHORT_PS_C2_T2","ONTSHORT_PS_C2_T3", "ONTSHORT_PS_C3_T1", "ONTSHORT_PS_C3_T2", "ONTSHORT_PS_C3_T3", "ONTSHORT_PS_T1_T2", "ONTSHORT_PS_T1_T3","ONTSHORT_PS_T2_T3",
                             "ONTLONG_AF_C1_C2", "ONTLONG_AF_C1_C3", "ONTLONG_AF_C1_T1", "ONTLONG_AF_C1_T2", "ONTLONG_AF_C1_T3", "ONTLONG_AF_C2_C3", "ONTLONG_AF_C2_T1", "ONTLONG_AF_C2_T2","ONTLONG_AF_C2_T3", "ONTLONG_AF_C3_T1", "ONTLONG_AF_C3_T2", "ONTLONG_AF_C3_T3", "ONTLONG_AF_T1_T2", "ONTLONG_AF_T1_T3","ONTLONG_AF_T2_T3",
                             "ONTLONG_PS_C1_C2", "ONTLONG_PS_C1_C3", "ONTLONG_PS_C1_T1", "ONTLONG_PS_C1_T2", "ONTLONG_PS_C1_T3", "ONTLONG_PS_C2_C3", "ONTLONG_PS_C2_T1", "ONTLONG_PS_C2_T2","ONTLONG_PS_C2_T3", "ONTLONG_PS_C3_T1", "ONTLONG_PS_C3_T2", "ONTLONG_PS_C3_T3", "ONTLONG_PS_T1_T2", "ONTLONG_PS_T1_T3","ONTLONG_PS_T2_T3")

beta_18s_bray_same_sampling = beta_18s_bray_same_sampling %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% mutate(method = str_extract(sample, "ILLUMINA|ONTSHORT|ONTLONG"), sampling = str_extract(sample, "AF|PS"), method_sampling = paste(method, sampling, sep = "_"), between = str_sub(sample, start = -8))
beta_18s_bray_same_sampling = beta_18s_bray_same_sampling %>% rowwise() %>% mutate(mean = mean(c_across(starts_with("value"))), sd = sd(c_across(starts_with("value"))), sem  = sd / sqrt(9)) %>% ungroup() 

beta_18s_bray_same_sampling$method = factor (beta_18s_bray_same_sampling$method, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
beta_18s_bray_same_sampling$sampling = factor (beta_18s_bray_same_sampling$sampling)
beta_18s_bray_same_sampling$method_sampling = factor (beta_18s_bray_same_sampling$method_sampling)
beta_18s_bray_same_sampling$between = factor (beta_18s_bray_same_sampling$between)
str(beta_18s_bray_same_sampling)

beta_18s_bray_same_env = beta_18s_bray_same_sampling %>% filter(str_sub(between, 4, 4) == str_sub(between, -2, -2))
beta_18s_bray_same_env = beta_18s_bray_same_env %>% mutate(environment_sampling = str_sub(between, 2, 4))
beta_18s_bray_same_env = beta_18s_bray_same_env %>% mutate(environment = str_sub(between, 4, 4))
beta_18s_bray_same_env = beta_18s_bray_same_env %>%
  pivot_longer(cols = starts_with("value"), 
               names_to = "value_id", 
               values_to = "value") %>%
  select(sample, method, sampling, method_sampling, between, mean, sd, sem, environment_sampling, environment, value_id, value)

beta_18s_bray_same_env$environment = factor(beta_18s_bray_same_env$environment, levels = c("C","T"))
str(beta_18s_bray_same_env)


# stats
shapiro.test(beta_16s_jaccard_same_env_average$mean)
shapiro.test(beta_16s_bray_same_env_average$mean) #

shapiro.test(beta_18s_jaccard_same_env_average$mean) 
shapiro.test(beta_18s_bray_same_env_average$mean) 

library(car)
leveneTest(mean ~ method, data = beta_16s_jaccard_same_env_average) 
leveneTest(mean ~ method, data = beta_16s_bray_same_env_average) 

leveneTest(mean ~ method, data = beta_18s_jaccard_same_env_average) 
leveneTest(mean ~ method, data = beta_18s_bray_same_env_average)

library(dunn.test)
kruskal.test(mean ~ method, data = beta_16s_jaccard_same_env_average)
dunn.test(beta_16s_jaccard_same_env_average$mean, beta_16s_jaccard_same_env_average$method, method = "bonferroni")


kruskal.test(mean ~ method, data = beta_16s_bray_same_env_average)
dunn.test(beta_16s_bray_same_env_average$mean, beta_16s_bray_same_env_average$method, method = "bonferroni")


kruskal.test(mean ~ method, data = beta_18s_jaccard_same_env_average)
dunn.test(beta_18s_jaccard_same_env_average$mean, beta_18s_jaccard_same_env_average$method, method = "bonferroni")


kruskal.test(mean ~ method, data = beta_18s_bray_same_env_average) 
dunn.test(beta_18s_bray_same_env_average$mean, beta_18s_bray_same_env_average$method, method = "bonferroni")


summary(beta_16s_jaccard_same_env_average[beta_16s_jaccard_same_env_average$method == "ILLUMINA", ]$mean)
summary(beta_16s_jaccard_same_env_average[beta_16s_jaccard_same_env_average$method == "ONTSHORT", ]$mean)
summary(beta_16s_jaccard_same_env_average[beta_16s_jaccard_same_env_average$method == "ONTLONG", ]$mean)

summary(beta_16s_bray_same_env_average[beta_16s_bray_same_env_average$method == "ILLUMINA", ]$mean)
summary(beta_16s_bray_same_env_average[beta_16s_bray_same_env_average$method == "ONTSHORT", ]$mean)
summary(beta_16s_bray_same_env_average[beta_16s_bray_same_env_average$method == "ONTLONG", ]$mean)

summary(beta_18s_jaccard_same_env_average[beta_18s_jaccard_same_env_average$method == "ILLUMINA", ]$mean)
summary(beta_18s_jaccard_same_env_average[beta_18s_jaccard_same_env_average$method == "ONTSHORT", ]$mean)
summary(beta_18s_jaccard_same_env_average[beta_18s_jaccard_same_env_average$method == "ONTLONG", ]$mean)

summary(beta_18s_bray_same_env_average[beta_18s_bray_same_env_average$method == "ILLUMINA", ]$mean)
summary(beta_18s_bray_same_env_average[beta_18s_bray_same_env_average$method == "ONTSHORT", ]$mean)
summary(beta_18s_bray_same_env_average[beta_18s_bray_same_env_average$method == "ONTLONG", ]$mean)


#  plots
library(tidyverse)
library(ggsignif)

beta_panel_16s_jaccard_average = ggplot(beta_16s_jaccard_same_env_average, aes(x = method, y = mean)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  ggtitle("") + xlab("") + ylab("Jaccard") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

beta_panel_16s_bray_average = ggplot(beta_16s_bray_same_env_average, aes(x = method, y = mean)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  ggtitle("") + xlab("") + ylab("Bray-Curtis") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

beta_panel_18s_jaccard_average = ggplot(beta_18s_jaccard_same_env_average, aes(x = method, y = mean)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  ggtitle("") + xlab("") + ylab("Jaccard") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

beta_panel_18s_bray_average = ggplot(beta_18s_bray_same_env_average, aes(x = method, y = mean)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(environment_sampling), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  ggtitle("") + xlab("") + ylab("Bray-Curtis") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 



















### gamma 
library(tidyverse)
library(vegan)

#data
otu_illu_16sv4_combined = otu_illu_16sv4_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")

otu_ont_16sv4_combined = otu_ont_16sv4_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")

otu_ont_fl16s_combined = otu_ont_fl16s_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")

otu_illu_18sv9_combined = otu_illu_18sv9_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")

otu_ont_18sv9_combined = otu_ont_18sv9_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")

otu_ont_fl18s_combined = otu_ont_fl18s_rarefied %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T",)) %>%
  group_by(group) %>% summarise(across(starts_with("otu"), sum), .groups = "drop") %>% column_to_rownames(var = "group")


gamma_OTU_16s = data.frame(sample = c(rownames(otu_illu_16sv4_combined), rownames(otu_ont_16sv4_combined), rownames(otu_ont_fl16s_combined)), 
                           value = c(specnumber(otu_illu_16sv4_combined), specnumber(otu_ont_16sv4_combined), specnumber(otu_ont_fl16s_combined)), 
                           method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 4))
gamma_OTU_16s$method = factor (gamma_OTU_16s$method, levels = c("Illumina","ONT_short","ONT_long"))

str(gamma_OTU_16s)


gamma_shannon_16s = data.frame(sample = c(rownames(otu_illu_16sv4_combined), rownames(otu_ont_16sv4_combined), rownames(otu_ont_fl16s_combined)), 
                              value = c(diversity(otu_illu_16sv4_combined, index = "shannon"), diversity(otu_ont_16sv4_combined, index = "shannon"), diversity(otu_ont_fl16s_combined, index = "shannon")), 
                              method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 4))
gamma_shannon_16s$method = factor (gamma_shannon_16s$method, levels = c("Illumina","ONT_short","ONT_long"))
str(gamma_shannon_16s)



gamma_OTU_18s = data.frame(sample = c(rownames(otu_illu_18sv9_combined), rownames(otu_ont_18sv9_combined), rownames(otu_ont_fl18s_combined)), 
                           value = c(specnumber(otu_illu_18sv9_combined), specnumber(otu_ont_18sv9_combined), specnumber(otu_ont_fl18s_combined)), 
                           method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 4))
gamma_OTU_18s$method = factor (gamma_OTU_18s$method, levels = c("Illumina","ONT_short","ONT_long"))
str(gamma_OTU_18s)


gamma_shannon_18s = data.frame(sample = c(rownames(otu_illu_18sv9_combined), rownames(otu_ont_18sv9_combined), rownames(otu_ont_fl18s_combined)), 
                               value = c(diversity(otu_illu_18sv9_combined, index = "shannon"), diversity(otu_ont_18sv9_combined, index = "shannon"), diversity(otu_ont_fl18s_combined, index = "shannon")), 
                               method = rep(c("Illumina", "ONT_short", "ONT_long"), each = 4))
gamma_shannon_18s$method = factor (gamma_shannon_18s$method, levels = c("Illumina","ONT_short","ONT_long"))
str(gamma_shannon_18s)


## stats
shapiro.test(gamma_OTU_16s$value) 
shapiro.test(gamma_shannon_16s$value) 

shapiro.test(gamma_OTU_18s$value)
shapiro.test(gamma_shannon_18s$value) 

library(car)
leveneTest(value ~ method, data = gamma_OTU_16s) 
leveneTest(value ~ method, data = gamma_shannon_16s) 

leveneTest(value ~ method, data = gamma_OTU_18s)
leveneTest(value ~ method, data = gamma_shannon_18s) 


library(dunn.test)
kruskal.test(value ~ method, data = gamma_OTU_16s) 
dunn.test(gamma_OTU_16s$value, gamma_OTU_16s$method, method = "bonferroni")

kruskal.test(value ~ method, data = gamma_shannon_16s) 
dunn.test(gamma_shannon_16s$value, gamma_shannon_16s$method, method = "bonferroni")


kruskal.test(value ~ method, data = gamma_OTU_18s) 
dunn.test(gamma_OTU_18s$value, gamma_OTU_18s$method, method = "bonferroni")

kruskal.test(value ~ method, data = gamma_shannon_18s)
dunn.test(gamma_shannon_18s$value, gamma_shannon_18s$method, method = "bonferroni")


##  plots
library(tidyverse)
library(ggsignif)

gamma_OTU_16s_figure = ggplot(gamma_OTU_16s, aes(x = method, y = value)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(sample), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 5000), breaks = c(0, 2500, 5000)) +
  ggtitle("") + xlab("") + ylab("Richness") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

gamma_OTU_18s_figure = ggplot(gamma_OTU_18s, aes(x = method, y = value)) +
  geom_point(size = 2, stroke = 0.5, aes(shape = factor(sample), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) +   
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) + 
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(0, 1000), breaks = c(0, 500, 1000)) +
  ggtitle("") + xlab("") + ylab("Richness") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

gamma_shannon_16s_figure = ggplot(gamma_shannon_16s, aes(x = method, y = value)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(sample), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(1, 7), breaks = c(1, 4, 7)) +
  ggtitle("") + xlab("") + ylab("Shannon") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

gamma_shannon_18s_figure = ggplot(gamma_shannon_18s, aes(x = method, y = value)) +
  geom_point(size = 2,  stroke = 0.5, aes(shape = factor(sample), color = factor(method)), position = position_dodge(width = 1), show.legend = F, alpha = 0.8) +
  scale_shape_manual(values = c("F_C" = 1, "S_C" = 0, "F_T" = 16, "S_T" = 15)) + 
  scale_color_manual(values = c("Illumina" = "#A9A9A9", "ONT_short" = "#9BC4E4", "ONT_long" = "#3579A1")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted", color = "black", linewidth = 0.3) +
  scale_x_discrete(label = c("Illumina","ONT Short", "ONT Long")) +
  scale_y_continuous(limits = c(1, 7), breaks = c(1, 4, 7)) +
  ggtitle("") + xlab("") + ylab("Shannon") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 




##  figure
library(patchwork)
figure_diversity = alpha_OTU_16s_average_figure + alpha_shannon_16s_average_figure + alpha_OTU_18s_average_figure + alpha_shannon_18s_average_figure + 
                   beta_panel_16s_jaccard_average + beta_panel_16s_bray_average + beta_panel_18s_jaccard_average + beta_panel_18s_bray_average + 
                   gamma_OTU_16s_figure + gamma_shannon_16s_figure  + gamma_OTU_18s_figure  + gamma_shannon_18s_figure +
                   plot_layout(ncol = 3, nrow = 4, byrow = F)

ggsave("figure_diversity.pdf", plot = figure_diversity, width = 17, height = 18, units = "cm", device = "pdf")














































### COMMUNITY COMPOSITION --------------------------------------------------------------------------------
library(tidyverse)

## data
# add taxon names
cc_illu_16sv4 = as.data.frame(otu_illu_16sv4_rarefied)
colnames(cc_illu_16sv4) = sapply(colnames(cc_illu_16sv4), function(x) {
  if (x %in% names(setNames(tax_illu_16sv4$Taxon,tax_illu_16sv4$OTU.ID))) {
    return(setNames(tax_illu_16sv4$Taxon,tax_illu_16sv4$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_illu_16sv4),100)

cc_ont_16sv4 = as.data.frame(otu_ont_16sv4_rarefied)
colnames(cc_ont_16sv4) = sapply(colnames(cc_ont_16sv4), function(x) {
  if (x %in% names(setNames(tax_ont_16sv4$Taxon,tax_ont_16sv4$OTU.ID))) {
    return(setNames(tax_ont_16sv4$Taxon,tax_ont_16sv4$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_ont_16sv4),100)

cc_ont_fl16s = as.data.frame(otu_ont_fl16s_rarefied)
colnames(cc_ont_fl16s) = sapply(colnames(cc_ont_fl16s), function(x) {
  if (x %in% names(setNames(tax_ont_fl16s$Taxon,tax_ont_fl16s$OTU.ID))) {
    return(setNames(tax_ont_fl16s$Taxon,tax_ont_fl16s$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_ont_fl16s),100)

cc_illu_18sv9 = as.data.frame(otu_illu_18sv9_rarefied)
colnames(cc_illu_18sv9) = sapply(colnames(cc_illu_18sv9), function(x) {
  if (x %in% names(setNames(tax_illu_18sv9$Taxon,tax_illu_18sv9$OTU.ID))) {
    return(setNames(tax_illu_18sv9$Taxon,tax_illu_18sv9$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_illu_18sv9),100)

cc_ont_18sv9 = as.data.frame(otu_ont_18sv9_rarefied)
colnames(cc_ont_18sv9) = sapply(colnames(cc_ont_18sv9), function(x) {
  if (x %in% names(setNames(tax_ont_18sv9$Taxon,tax_ont_18sv9$OTU.ID))) {
    return(setNames(tax_ont_18sv9$Taxon,tax_ont_18sv9$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_ont_18sv9),100)

cc_ont_fl18s = as.data.frame(otu_ont_fl18s_rarefied)
colnames(cc_ont_fl18s) = sapply(colnames(cc_ont_fl18s), function(x) {
  if (x %in% names(setNames(tax_ont_fl18s$Taxon,tax_ont_fl18s$OTU.ID))) {
    return(setNames(tax_ont_fl18s$Taxon,tax_ont_fl18s$OTU.ID)[x])
  } else {
    return(x) 
  }
})
head(colnames(cc_ont_fl18s),100)


# define extract functions
extract_phylum_name = function(x) {
  match = str_extract(x, "(?<=p__)[^;]+") 
  return(match)
}

extract_dvs_name = function(x) {
  match = str_extract(x, "(?<=dvs__)[^;]+") 
  return(match)
}

extract_genus_name = function(x) {
  match = str_extract(x, "(?<=g__)[^;]+") 
  return(match)
}


# top phylum 16s
rownames(cc_illu_16sv4_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())
rownames(cc_ont_16sv4_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())
rownames(cc_ont_fl16s_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())

intersect(intersect(rownames(cc_illu_16sv4_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame()), rownames(cc_ont_16sv4_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())), rownames(cc_ont_fl16s_phylum_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame()))

# top subdivision 18s
rownames(cc_illu_18sv9_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())
rownames(cc_ont_18sv9_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())
rownames(cc_ont_fl18s_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())

intersect(intersect(rownames(cc_illu_18sv9_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame()), rownames(cc_ont_18sv9_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame())), rownames(cc_ont_fl18s_subdivision_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:20] %>% data.frame()))

# top genus 16s
rownames(cc_illu_16sv4_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())
rownames(cc_ont_16sv4_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())
rownames(cc_ont_fl16s_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())

intersect(intersect(rownames(cc_illu_16sv4_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame()), rownames(cc_ont_16sv4_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())), rownames(cc_ont_fl16s_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame()))

# top genus 18s
rownames(cc_illu_18sv9_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())
rownames(cc_ont_18sv9_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())
rownames(cc_ont_fl18s_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())

intersect(intersect(rownames(cc_illu_18sv9_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame()), rownames(cc_ont_18sv9_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame())), rownames(cc_ont_fl18s_genus_cleaned %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:30] %>% data.frame()))




cc_illu_16sv4_phylum = cc_illu_16sv4
colnames(cc_illu_16sv4_phylum) = sapply(colnames(cc_illu_16sv4_phylum), extract_phylum_name)
head(colnames(cc_illu_16sv4_phylum),100)
sum(is.na(colnames(cc_illu_16sv4_phylum))) 

cc_illu_16sv4_phylum = cc_illu_16sv4_phylum[ , !is.na(colnames(cc_illu_16sv4_phylum))]
head(colnames(cc_illu_16sv4_phylum),100)

colnames(cc_illu_16sv4_phylum) = make.unique(colnames(cc_illu_16sv4_phylum))
head(colnames(cc_illu_16sv4_phylum),100)

cc_illu_16sv4_phylum = as.data.frame(cc_illu_16sv4_phylum / rowSums(cc_illu_16sv4_phylum) *100)
rowSums(cc_illu_16sv4_phylum) 

cc_illu_16sv4_phylum_cleaned = cc_illu_16sv4_phylum %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_illu_16sv4_phylum_cleaned = cc_illu_16sv4_phylum_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_illu_16sv4_phylum_cleaned = cc_illu_16sv4_phylum_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_illu_16sv4_phylum_cleaned = cc_illu_16sv4_phylum_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_illu_16sv4_phylum_cleaned = as.data.frame(cc_illu_16sv4_phylum_cleaned)
rownames(cc_illu_16sv4_phylum_cleaned) = cc_illu_16sv4_phylum_cleaned$Sample
cc_illu_16sv4_phylum_cleaned = cc_illu_16sv4_phylum_cleaned[, -1] 
rowSums(cc_illu_16sv4_phylum_cleaned) 
colnames(cc_illu_16sv4_phylum_cleaned) 

cc_illu_16sv4_phylum_assigned = cc_illu_16sv4_phylum_cleaned %>% select("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                        "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                        "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543")
sum(cc_illu_16sv4_phylum_assigned) 
colnames(cc_illu_16sv4_phylum_assigned) 

cc_illu_16sv4_phylum_assigned = cc_illu_16sv4_phylum_assigned %>%
  mutate("Others" = 100 - rowSums(cc_illu_16sv4_phylum_assigned))

cc_illu_16sv4_phylum_assigned_avg = cc_illu_16sv4_phylum_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_illu_16sv4_phylum_assigned_avg) 

cc_illu_16sv4_phylum_assigned_avg = cc_illu_16sv4_phylum_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_illu_16sv4_phylum_assigned_avg[,"Proportion"]) 

cc_illu_16sv4_phylum_assigned_avg$Sample = factor(cc_illu_16sv4_phylum_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                     "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                     "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                     "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_illu_16sv4_phylum_assigned_avg$Group = factor(cc_illu_16sv4_phylum_assigned_avg$Group, levels = c("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                                                     "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                                                     "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543",
                                                                                                     "Others"))
str(cc_illu_16sv4_phylum_assigned_avg)




cc_ont_16sv4_phylum = cc_ont_16sv4
colnames(cc_ont_16sv4_phylum) = sapply(colnames(cc_ont_16sv4_phylum), extract_phylum_name)
head(colnames(cc_ont_16sv4_phylum),100)
sum(is.na(colnames(cc_ont_16sv4_phylum))) 

cc_ont_16sv4_phylum = cc_ont_16sv4_phylum[ , !is.na(colnames(cc_ont_16sv4_phylum))]
head(colnames(cc_ont_16sv4_phylum),100)

colnames(cc_ont_16sv4_phylum) = make.unique(colnames(cc_ont_16sv4_phylum))
head(colnames(cc_ont_16sv4_phylum),100)

cc_ont_16sv4_phylum = as.data.frame(cc_ont_16sv4_phylum / rowSums(cc_ont_16sv4_phylum) *100)
rowSums(cc_ont_16sv4_phylum)

cc_ont_16sv4_phylum_cleaned = cc_ont_16sv4_phylum %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_16sv4_phylum_cleaned = cc_ont_16sv4_phylum_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_16sv4_phylum_cleaned = cc_ont_16sv4_phylum_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_16sv4_phylum_cleaned = cc_ont_16sv4_phylum_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_16sv4_phylum_cleaned = as.data.frame(cc_ont_16sv4_phylum_cleaned)
rownames(cc_ont_16sv4_phylum_cleaned) = cc_ont_16sv4_phylum_cleaned$Sample
cc_ont_16sv4_phylum_cleaned = cc_ont_16sv4_phylum_cleaned[, -1] 
rowSums(cc_ont_16sv4_phylum_cleaned) 
colnames(cc_ont_16sv4_phylum_cleaned) 

cc_ont_16sv4_phylum_assigned = cc_ont_16sv4_phylum_cleaned %>% select("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                        "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                        "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543")
sum(cc_ont_16sv4_phylum_assigned) 
colnames(cc_ont_16sv4_phylum_assigned) 

cc_ont_16sv4_phylum_assigned = cc_ont_16sv4_phylum_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_16sv4_phylum_assigned))

cc_ont_16sv4_phylum_assigned_avg = cc_ont_16sv4_phylum_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_16sv4_phylum_assigned_avg) 

cc_ont_16sv4_phylum_assigned_avg = cc_ont_16sv4_phylum_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_16sv4_phylum_assigned_avg[,"Proportion"]) 

cc_ont_16sv4_phylum_assigned_avg$Sample = factor(cc_ont_16sv4_phylum_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                     "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                     "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                     "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_16sv4_phylum_assigned_avg$Group = factor(cc_ont_16sv4_phylum_assigned_avg$Group, levels = c("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                                                     "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                                                     "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543",
                                                                                                     "Others"))
str(cc_ont_16sv4_phylum_assigned_avg)




cc_ont_fl16s_phylum = cc_ont_fl16s
colnames(cc_ont_fl16s_phylum) = sapply(colnames(cc_ont_fl16s_phylum), extract_phylum_name)
head(colnames(cc_ont_fl16s_phylum),100)
sum(is.na(colnames(cc_ont_fl16s_phylum)))

cc_ont_fl16s_phylum = cc_ont_fl16s_phylum[ , !is.na(colnames(cc_ont_fl16s_phylum))]
head(colnames(cc_ont_fl16s_phylum),100)

colnames(cc_ont_fl16s_phylum) = make.unique(colnames(cc_ont_fl16s_phylum))
head(colnames(cc_ont_fl16s_phylum),100)

cc_ont_fl16s_phylum = as.data.frame(cc_ont_fl16s_phylum / rowSums(cc_ont_fl16s_phylum) *100)
rowSums(cc_ont_fl16s_phylum)

cc_ont_fl16s_phylum_cleaned = cc_ont_fl16s_phylum %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_fl16s_phylum_cleaned = cc_ont_fl16s_phylum_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_fl16s_phylum_cleaned = cc_ont_fl16s_phylum_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_fl16s_phylum_cleaned = cc_ont_fl16s_phylum_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_fl16s_phylum_cleaned = as.data.frame(cc_ont_fl16s_phylum_cleaned)
rownames(cc_ont_fl16s_phylum_cleaned) = cc_ont_fl16s_phylum_cleaned$Sample
cc_ont_fl16s_phylum_cleaned = cc_ont_fl16s_phylum_cleaned[, -1] 
rowSums(cc_ont_fl16s_phylum_cleaned) 
colnames(cc_ont_fl16s_phylum_cleaned) 

cc_ont_fl16s_phylum_assigned = cc_ont_fl16s_phylum_cleaned %>% select("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                        "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                        "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543")
sum(cc_ont_fl16s_phylum_assigned) 
colnames(cc_ont_fl16s_phylum_assigned) 

cc_ont_fl16s_phylum_assigned = cc_ont_fl16s_phylum_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_fl16s_phylum_assigned))

cc_ont_fl16s_phylum_assigned_avg = cc_ont_fl16s_phylum_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_fl16s_phylum_assigned_avg) 

cc_ont_fl16s_phylum_assigned_avg = cc_ont_fl16s_phylum_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_fl16s_phylum_assigned_avg[,"Proportion"]) 

cc_ont_fl16s_phylum_assigned_avg$Sample = factor(cc_ont_fl16s_phylum_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                       "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                       "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                       "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_fl16s_phylum_assigned_avg$Group = factor(cc_ont_fl16s_phylum_assigned_avg$Group, levels = c("Pseudomonadota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Planctomycetota",
                                                                                                     "Verrucomicrobiota", "Marinisomatota", "Myxococcota_A_473307", "SAR324", "Bacillota_I",
                                                                                                     "Fusobacteriota", "Bdellovibrionota_473306", "Campylobacterota_A", "Thermosulfidibacterota", "Desulfobacterota_G_459543",
                                                                                                     "Others"))
str(cc_ont_fl16s_phylum_assigned_avg)







cc_illu_18sv9_subdivision = cc_illu_18sv9
colnames(cc_illu_18sv9_subdivision) = sapply(colnames(cc_illu_18sv9_subdivision), extract_dvs_name)
head(colnames(cc_illu_18sv9_subdivision),100)
sum(is.na(colnames(cc_illu_18sv9_subdivision))) 

cc_illu_18sv9_subdivision = cc_illu_18sv9_subdivision[ , !is.na(colnames(cc_illu_18sv9_subdivision))]
head(colnames(cc_illu_18sv9_subdivision),100)

colnames(cc_illu_18sv9_subdivision) = make.unique(colnames(cc_illu_18sv9_subdivision))
head(colnames(cc_illu_18sv9_subdivision),100)

cc_illu_18sv9_subdivision = as.data.frame(cc_illu_18sv9_subdivision / rowSums(cc_illu_18sv9_subdivision) *100)
rowSums(cc_illu_18sv9_subdivision) 

cc_illu_18sv9_subdivision_cleaned = cc_illu_18sv9_subdivision %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_illu_18sv9_subdivision_cleaned = cc_illu_18sv9_subdivision_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_illu_18sv9_subdivision_cleaned = cc_illu_18sv9_subdivision_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_illu_18sv9_subdivision_cleaned = cc_illu_18sv9_subdivision_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_illu_18sv9_subdivision_cleaned = as.data.frame(cc_illu_18sv9_subdivision_cleaned)
rownames(cc_illu_18sv9_subdivision_cleaned) = cc_illu_18sv9_subdivision_cleaned$Sample
cc_illu_18sv9_subdivision_cleaned = cc_illu_18sv9_subdivision_cleaned[, -1] 
rowSums(cc_illu_18sv9_subdivision_cleaned) 
colnames(cc_illu_18sv9_subdivision_cleaned) 

cc_illu_18sv9_subdivision_assigned = cc_illu_18sv9_subdivision_cleaned %>% select("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                  "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                  "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida")
sum(cc_illu_18sv9_subdivision_assigned) 
colnames(cc_illu_18sv9_subdivision_assigned) 

cc_illu_18sv9_subdivision_assigned = cc_illu_18sv9_subdivision_assigned %>%
  mutate("Others" = 100 - rowSums(cc_illu_18sv9_subdivision_assigned))

cc_illu_18sv9_subdivision_assigned_avg = cc_illu_18sv9_subdivision_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_illu_18sv9_subdivision_assigned_avg) 

cc_illu_18sv9_subdivision_assigned_avg = cc_illu_18sv9_subdivision_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_illu_18sv9_subdivision_assigned_avg[,"Proportion"]) 

cc_illu_18sv9_subdivision_assigned_avg$Sample = factor(cc_illu_18sv9_subdivision_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                       "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                       "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                       "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_illu_18sv9_subdivision_assigned_avg$Group = factor(cc_illu_18sv9_subdivision_assigned_avg$Group, levels = c("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                                               "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                                               "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida",
                                                                                                               "Others"))
str(cc_illu_18sv9_subdivision_assigned_avg)







cc_ont_18sv9_subdivision = cc_ont_18sv9
colnames(cc_ont_18sv9_subdivision) = sapply(colnames(cc_ont_18sv9_subdivision), extract_dvs_name)
head(colnames(cc_ont_18sv9_subdivision),100)
sum(is.na(colnames(cc_ont_18sv9_subdivision))) 

cc_ont_18sv9_subdivision = cc_ont_18sv9_subdivision[ , !is.na(colnames(cc_ont_18sv9_subdivision))]
head(colnames(cc_ont_18sv9_subdivision),100)

colnames(cc_ont_18sv9_subdivision) = make.unique(colnames(cc_ont_18sv9_subdivision))
head(colnames(cc_ont_18sv9_subdivision),100)

cc_ont_18sv9_subdivision = as.data.frame(cc_ont_18sv9_subdivision / rowSums(cc_ont_18sv9_subdivision) *100)
rowSums(cc_ont_18sv9_subdivision) 

cc_ont_18sv9_subdivision_cleaned = cc_ont_18sv9_subdivision %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_18sv9_subdivision_cleaned = cc_ont_18sv9_subdivision_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_18sv9_subdivision_cleaned = cc_ont_18sv9_subdivision_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_18sv9_subdivision_cleaned = cc_ont_18sv9_subdivision_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_18sv9_subdivision_cleaned = as.data.frame(cc_ont_18sv9_subdivision_cleaned)
rownames(cc_ont_18sv9_subdivision_cleaned) = cc_ont_18sv9_subdivision_cleaned$Sample
cc_ont_18sv9_subdivision_cleaned = cc_ont_18sv9_subdivision_cleaned[, -1] 
rowSums(cc_ont_18sv9_subdivision_cleaned) 
colnames(cc_ont_18sv9_subdivision_cleaned) 

cc_ont_18sv9_subdivision_assigned = cc_ont_18sv9_subdivision_cleaned %>% select("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                  "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                  "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida")
sum(cc_ont_18sv9_subdivision_assigned) 
colnames(cc_ont_18sv9_subdivision_assigned) 

cc_ont_18sv9_subdivision_assigned = cc_ont_18sv9_subdivision_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_18sv9_subdivision_assigned))

cc_ont_18sv9_subdivision_assigned_avg = cc_ont_18sv9_subdivision_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_18sv9_subdivision_assigned_avg) 

cc_ont_18sv9_subdivision_assigned_avg = cc_ont_18sv9_subdivision_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_18sv9_subdivision_assigned_avg[,"Proportion"]) 

cc_ont_18sv9_subdivision_assigned_avg$Sample = factor(cc_ont_18sv9_subdivision_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                                 "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                                 "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                                 "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_18sv9_subdivision_assigned_avg$Group = factor(cc_ont_18sv9_subdivision_assigned_avg$Group, levels = c("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                                               "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                                               "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida",
                                                                                                               "Others"))
str(cc_ont_18sv9_subdivision_assigned_avg)



cc_ont_fl18s_subdivision = cc_ont_fl18s
colnames(cc_ont_fl18s_subdivision) = sapply(colnames(cc_ont_fl18s_subdivision), extract_dvs_name)
head(colnames(cc_ont_fl18s_subdivision),100)
sum(is.na(colnames(cc_ont_fl18s_subdivision))) 

cc_ont_fl18s_subdivision = cc_ont_fl18s_subdivision[ , !is.na(colnames(cc_ont_fl18s_subdivision))]
head(colnames(cc_ont_fl18s_subdivision),100)

colnames(cc_ont_fl18s_subdivision) = make.unique(colnames(cc_ont_fl18s_subdivision))
head(colnames(cc_ont_fl18s_subdivision),100)

cc_ont_fl18s_subdivision = as.data.frame(cc_ont_fl18s_subdivision / rowSums(cc_ont_fl18s_subdivision) *100)
rowSums(cc_ont_fl18s_subdivision)

cc_ont_fl18s_subdivision_cleaned = cc_ont_fl18s_subdivision %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_fl18s_subdivision_cleaned = cc_ont_fl18s_subdivision_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_fl18s_subdivision_cleaned = cc_ont_fl18s_subdivision_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_fl18s_subdivision_cleaned = cc_ont_fl18s_subdivision_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_fl18s_subdivision_cleaned = as.data.frame(cc_ont_fl18s_subdivision_cleaned)
rownames(cc_ont_fl18s_subdivision_cleaned) = cc_ont_fl18s_subdivision_cleaned$Sample
cc_ont_fl18s_subdivision_cleaned = cc_ont_fl18s_subdivision_cleaned[, -1] 
rowSums(cc_ont_fl18s_subdivision_cleaned) 
colnames(cc_ont_fl18s_subdivision_cleaned) 

cc_ont_fl18s_subdivision_assigned = cc_ont_fl18s_subdivision_cleaned %>% select("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida")
sum(cc_ont_fl18s_subdivision_assigned) 
colnames(cc_ont_fl18s_subdivision_assigned) 

cc_ont_fl18s_subdivision_assigned = cc_ont_fl18s_subdivision_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_fl18s_subdivision_assigned))

cc_ont_fl18s_subdivision_assigned_avg = cc_ont_fl18s_subdivision_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_fl18s_subdivision_assigned_avg) 

cc_ont_fl18s_subdivision_assigned_avg = cc_ont_fl18s_subdivision_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_fl18s_subdivision_assigned_avg[,"Proportion"]) 

cc_ont_fl18s_subdivision_assigned_avg$Sample = factor(cc_ont_fl18s_subdivision_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                               "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                               "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                               "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_fl18s_subdivision_assigned_avg$Group = factor(cc_ont_fl18s_subdivision_assigned_avg$Group, levels = c("Gyrista", "Dinoflagellata", "Ciliophora", "Bigyra", "Chlorophyta_X",
                                                                                                             "Cryptophyta_X", "Cercozoa", "Haptophyta_X", "Picozoa_X", "Rhodophyta_X",
                                                                                                             "Radiolaria", "Fungi", "Ancyromonadida_X", "Apusomonada_X", "Kathablepharida",
                                                                                                             "Others"))
str(cc_ont_fl18s_subdivision_assigned_avg)








cc_illu_16sv4_genus = cc_illu_16sv4
colnames(cc_illu_16sv4_genus) = sapply(colnames(cc_illu_16sv4_genus), extract_genus_name)
head(colnames(cc_illu_16sv4_genus),100)
sum(is.na(colnames(cc_illu_16sv4_genus))) 

cc_illu_16sv4_genus = cc_illu_16sv4_genus[ , !is.na(colnames(cc_illu_16sv4_genus))]
head(colnames(cc_illu_16sv4_genus),100)

colnames(cc_illu_16sv4_genus) = make.unique(colnames(cc_illu_16sv4_genus))
head(colnames(cc_illu_16sv4_genus),100)

cc_illu_16sv4_genus = as.data.frame(cc_illu_16sv4_genus / rowSums(cc_illu_16sv4_genus) *100)
rowSums(cc_illu_16sv4_genus) 

cc_illu_16sv4_genus_cleaned = cc_illu_16sv4_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_illu_16sv4_genus_cleaned = cc_illu_16sv4_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_illu_16sv4_genus_cleaned = cc_illu_16sv4_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_illu_16sv4_genus_cleaned = cc_illu_16sv4_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_illu_16sv4_genus_cleaned = as.data.frame(cc_illu_16sv4_genus_cleaned)
rownames(cc_illu_16sv4_genus_cleaned) = cc_illu_16sv4_genus_cleaned$Sample
cc_illu_16sv4_genus_cleaned = cc_illu_16sv4_genus_cleaned[, -1] 
rowSums(cc_illu_16sv4_genus_cleaned)
colnames(cc_illu_16sv4_genus_cleaned) 

cc_illu_16sv4_genus_assigned = cc_illu_16sv4_genus_cleaned %>% select("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                      "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                      "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62")

sum(cc_illu_16sv4_genus_assigned) 
colnames(cc_illu_16sv4_genus_assigned) 

cc_illu_16sv4_genus_assigned = cc_illu_16sv4_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_illu_16sv4_genus_assigned))

cc_illu_16sv4_genus_assigned_avg = cc_illu_16sv4_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_illu_16sv4_genus_assigned_avg) 

cc_illu_16sv4_genus_assigned_avg = cc_illu_16sv4_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                                     "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                                     "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                                     "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_illu_16sv4_genus_assigned_avg[,"Proportion"]) 

cc_illu_16sv4_genus_assigned_avg$Sample = factor(cc_illu_16sv4_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                     "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                     "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                     "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_illu_16sv4_genus_assigned_avg$Group = factor(cc_illu_16sv4_genus_assigned_avg$Group, levels = c("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                                                   "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                                                   "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62",
                                                                                                   "Others"))
str(cc_illu_16sv4_genus_assigned_avg)




cc_ont_16sv4_genus = cc_ont_16sv4
colnames(cc_ont_16sv4_genus) = sapply(colnames(cc_ont_16sv4_genus), extract_genus_name)
head(colnames(cc_ont_16sv4_genus),100)
sum(is.na(colnames(cc_ont_16sv4_genus)))

cc_ont_16sv4_genus = cc_ont_16sv4_genus[ , !is.na(colnames(cc_ont_16sv4_genus))]
head(colnames(cc_ont_16sv4_genus),100)

colnames(cc_ont_16sv4_genus) = make.unique(colnames(cc_ont_16sv4_genus))
head(colnames(cc_ont_16sv4_genus),100)

cc_ont_16sv4_genus = as.data.frame(cc_ont_16sv4_genus / rowSums(cc_ont_16sv4_genus) *100)
rowSums(cc_ont_16sv4_genus) 

cc_ont_16sv4_genus_cleaned = cc_ont_16sv4_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_16sv4_genus_cleaned = cc_ont_16sv4_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_16sv4_genus_cleaned = cc_ont_16sv4_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_16sv4_genus_cleaned = cc_ont_16sv4_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_16sv4_genus_cleaned = as.data.frame(cc_ont_16sv4_genus_cleaned)
rownames(cc_ont_16sv4_genus_cleaned) = cc_ont_16sv4_genus_cleaned$Sample
cc_ont_16sv4_genus_cleaned = cc_ont_16sv4_genus_cleaned[, -1] 
rowSums(cc_ont_16sv4_genus_cleaned) 
colnames(cc_ont_16sv4_genus_cleaned) 

cc_ont_16sv4_genus_assigned = cc_ont_16sv4_genus_cleaned %>% select("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                    "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                    "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62")

sum(cc_ont_16sv4_genus_assigned)
colnames(cc_ont_16sv4_genus_assigned) 

cc_ont_16sv4_genus_assigned = cc_ont_16sv4_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_16sv4_genus_assigned))

cc_ont_16sv4_genus_assigned_avg = cc_ont_16sv4_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_16sv4_genus_assigned_avg) 

cc_ont_16sv4_genus_assigned_avg = cc_ont_16sv4_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                                     "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                                     "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                                     "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_16sv4_genus_assigned_avg[,"Proportion"]) 

cc_ont_16sv4_genus_assigned_avg$Sample = factor(cc_ont_16sv4_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                   "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                   "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                   "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_16sv4_genus_assigned_avg$Group = factor(cc_ont_16sv4_genus_assigned_avg$Group, levels = c("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                                                 "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                                                 "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62",
                                                                                                 "Others"))
str(cc_ont_16sv4_genus_assigned_avg)





cc_ont_fl16s_genus = cc_ont_fl16s
colnames(cc_ont_fl16s_genus) = sapply(colnames(cc_ont_fl16s_genus), extract_genus_name)
head(colnames(cc_ont_fl16s_genus),100)
sum(is.na(colnames(cc_ont_fl16s_genus))) 

cc_ont_fl16s_genus = cc_ont_fl16s_genus[ , !is.na(colnames(cc_ont_fl16s_genus))]
head(colnames(cc_ont_fl16s_genus),100) 

colnames(cc_ont_fl16s_genus) = make.unique(colnames(cc_ont_fl16s_genus))
head(colnames(cc_ont_fl16s_genus),100)

cc_ont_fl16s_genus = as.data.frame(cc_ont_fl16s_genus / rowSums(cc_ont_fl16s_genus) *100)
rowSums(cc_ont_fl16s_genus)

cc_ont_fl16s_genus_cleaned = cc_ont_fl16s_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_fl16s_genus_cleaned = cc_ont_fl16s_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_fl16s_genus_cleaned = cc_ont_fl16s_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_fl16s_genus_cleaned = cc_ont_fl16s_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_fl16s_genus_cleaned = as.data.frame(cc_ont_fl16s_genus_cleaned)
rownames(cc_ont_fl16s_genus_cleaned) = cc_ont_fl16s_genus_cleaned$Sample
cc_ont_fl16s_genus_cleaned = cc_ont_fl16s_genus_cleaned[, -1] 
rowSums(cc_ont_fl16s_genus_cleaned) 
colnames(cc_ont_fl16s_genus_cleaned) 

cc_ont_fl16s_genus_assigned = cc_ont_fl16s_genus_cleaned %>% select("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                    "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                    "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62")

sum(cc_ont_fl16s_genus_assigned)
colnames(cc_ont_fl16s_genus_assigned) 

cc_ont_fl16s_genus_assigned = cc_ont_fl16s_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_fl16s_genus_assigned))

cc_ont_fl16s_genus_assigned_avg = cc_ont_fl16s_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_fl16s_genus_assigned_avg) 

cc_ont_fl16s_genus_assigned_avg = cc_ont_fl16s_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                                     "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                                     "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                                     "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_fl16s_genus_assigned_avg[,"Proportion"]) 

cc_ont_fl16s_genus_assigned_avg$Sample = factor(cc_ont_fl16s_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                   "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                   "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                   "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_fl16s_genus_assigned_avg$Group = factor(cc_ont_fl16s_genus_assigned_avg$Group, levels = c("Marinagarivorans", "Actinomarina", "Parasynechococcus", "HIMB59", "Cellvibrio",
                                                                                                 "Luminiphilus", "MED-G82", "Pelagibacter_A_533952", "TMED14", "UBA10364",
                                                                                                 "UBA1275", "BACL14", "UBA1014", "UBA3031", "RS62",
                                                                                                 "Others"))
str(cc_ont_fl16s_genus_assigned_avg)





cc_illu_18sv9_genus = cc_illu_18sv9
colnames(cc_illu_18sv9_genus) = sapply(colnames(cc_illu_18sv9_genus), extract_genus_name)
head(colnames(cc_illu_18sv9_genus),100)
sum(is.na(colnames(cc_illu_18sv9_genus))) 

cc_illu_18sv9_genus = cc_illu_18sv9_genus[ , !is.na(colnames(cc_illu_18sv9_genus))]
head(colnames(cc_illu_18sv9_genus),100)

colnames(cc_illu_18sv9_genus) = make.unique(colnames(cc_illu_18sv9_genus))
head(colnames(cc_illu_18sv9_genus),100)

cc_illu_18sv9_genus = as.data.frame(cc_illu_18sv9_genus / rowSums(cc_illu_18sv9_genus) *100)
rowSums(cc_illu_18sv9_genus) 

cc_illu_18sv9_genus_cleaned = cc_illu_18sv9_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_illu_18sv9_genus_cleaned = cc_illu_18sv9_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_illu_18sv9_genus_cleaned = cc_illu_18sv9_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_illu_18sv9_genus_cleaned = cc_illu_18sv9_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_illu_18sv9_genus_cleaned = as.data.frame(cc_illu_18sv9_genus_cleaned)
rownames(cc_illu_18sv9_genus_cleaned) = cc_illu_18sv9_genus_cleaned$Sample
cc_illu_18sv9_genus_cleaned = cc_illu_18sv9_genus_cleaned[, -1] 
rowSums(cc_illu_18sv9_genus_cleaned) 
colnames(cc_illu_18sv9_genus_cleaned) 

cc_illu_18sv9_genus_assigned = cc_illu_18sv9_genus_cleaned %>% select("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                      "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                      "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia")

sum(cc_illu_18sv9_genus_assigned) 
colnames(cc_illu_18sv9_genus_assigned) 

cc_illu_18sv9_genus_assigned = cc_illu_18sv9_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_illu_18sv9_genus_assigned))

cc_illu_18sv9_genus_assigned_avg = cc_illu_18sv9_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_illu_18sv9_genus_assigned_avg) 

cc_illu_18sv9_genus_assigned_avg = cc_illu_18sv9_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_illu_18sv9_genus_assigned_avg[,"Proportion"]) # 1200

cc_illu_18sv9_genus_assigned_avg$Sample = factor(cc_illu_18sv9_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                     "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                     "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                     "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_illu_18sv9_genus_assigned_avg$Group = factor(cc_illu_18sv9_genus_assigned_avg$Group, levels = c("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                                                   "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                                                   "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia",
                                                                                                   "Others"))
str(cc_illu_18sv9_genus_assigned_avg)









cc_ont_18sv9_genus = cc_ont_18sv9
colnames(cc_ont_18sv9_genus) = sapply(colnames(cc_ont_18sv9_genus), extract_genus_name)
head(colnames(cc_ont_18sv9_genus),100)
sum(is.na(colnames(cc_ont_18sv9_genus))) 

cc_ont_18sv9_genus = cc_ont_18sv9_genus[ , !is.na(colnames(cc_ont_18sv9_genus))]
head(colnames(cc_ont_18sv9_genus),100)

colnames(cc_ont_18sv9_genus) = make.unique(colnames(cc_ont_18sv9_genus))
head(colnames(cc_ont_18sv9_genus),100)

cc_ont_18sv9_genus = as.data.frame(cc_ont_18sv9_genus / rowSums(cc_ont_18sv9_genus) *100)
rowSums(cc_ont_18sv9_genus) 

cc_ont_18sv9_genus_cleaned = cc_ont_18sv9_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_18sv9_genus_cleaned = cc_ont_18sv9_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_18sv9_genus_cleaned = cc_ont_18sv9_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_18sv9_genus_cleaned = cc_ont_18sv9_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_18sv9_genus_cleaned = as.data.frame(cc_ont_18sv9_genus_cleaned)
rownames(cc_ont_18sv9_genus_cleaned) = cc_ont_18sv9_genus_cleaned$Sample
cc_ont_18sv9_genus_cleaned = cc_ont_18sv9_genus_cleaned[, -1] 
rowSums(cc_ont_18sv9_genus_cleaned) 
colnames(cc_ont_18sv9_genus_cleaned) 

cc_ont_18sv9_genus_assigned = cc_ont_18sv9_genus_cleaned %>% select("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                      "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                      "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia")

sum(cc_ont_18sv9_genus_assigned) 
colnames(cc_ont_18sv9_genus_assigned) 

cc_ont_18sv9_genus_assigned = cc_ont_18sv9_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_18sv9_genus_assigned))

cc_ont_18sv9_genus_assigned_avg = cc_ont_18sv9_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_18sv9_genus_assigned_avg)

cc_ont_18sv9_genus_assigned_avg = cc_ont_18sv9_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_18sv9_genus_assigned_avg[,"Proportion"]) 

cc_ont_18sv9_genus_assigned_avg$Sample = factor(cc_ont_18sv9_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                     "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                     "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                     "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_18sv9_genus_assigned_avg$Group = factor(cc_ont_18sv9_genus_assigned_avg$Group, levels = c("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                                                   "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                                                   "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia",
                                                                                                   "Others"))
str(cc_ont_18sv9_genus_assigned_avg)





cc_ont_fl18s_genus = cc_ont_fl18s
colnames(cc_ont_fl18s_genus) = sapply(colnames(cc_ont_fl18s_genus), extract_genus_name)
head(colnames(cc_ont_fl18s_genus),100)
sum(is.na(colnames(cc_ont_fl18s_genus))) 

cc_ont_fl18s_genus = cc_ont_fl18s_genus[ , !is.na(colnames(cc_ont_fl18s_genus))]
head(colnames(cc_ont_fl18s_genus),100)

colnames(cc_ont_fl18s_genus) = make.unique(colnames(cc_ont_fl18s_genus))
head(colnames(cc_ont_fl18s_genus),100)

cc_ont_fl18s_genus = as.data.frame(cc_ont_fl18s_genus / rowSums(cc_ont_fl18s_genus) *100)
rowSums(cc_ont_fl18s_genus) 

cc_ont_fl18s_genus_cleaned = cc_ont_fl18s_genus %>% mutate(Sample = rownames(.)) %>% pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")
cc_ont_fl18s_genus_cleaned = cc_ont_fl18s_genus_cleaned %>% mutate(Group_unique = gsub("\\..*", "", Group))
cc_ont_fl18s_genus_cleaned = cc_ont_fl18s_genus_cleaned %>% group_by(Sample, Group_unique) %>% summarise(Value = sum(Value), .groups = "drop")
cc_ont_fl18s_genus_cleaned = cc_ont_fl18s_genus_cleaned %>% pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

cc_ont_fl18s_genus_cleaned = as.data.frame(cc_ont_fl18s_genus_cleaned)
rownames(cc_ont_fl18s_genus_cleaned) = cc_ont_fl18s_genus_cleaned$Sample
cc_ont_fl18s_genus_cleaned = cc_ont_fl18s_genus_cleaned[, -1] 
rowSums(cc_ont_fl18s_genus_cleaned) 
colnames(cc_ont_fl18s_genus_cleaned) 

cc_ont_fl18s_genus_assigned = cc_ont_fl18s_genus_cleaned %>% select("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                    "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                    "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia")

sum(cc_ont_fl18s_genus_assigned) 
colnames(cc_ont_fl18s_genus_assigned) 

cc_ont_fl18s_genus_assigned = cc_ont_fl18s_genus_assigned %>%
  mutate("Others" = 100 - rowSums(cc_ont_fl18s_genus_assigned))

cc_ont_fl18s_genus_assigned_avg = cc_ont_fl18s_genus_assigned %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = case_when(
    grepl("^F_C1_", sample) ~ "F_C1", grepl("^F_C2_", sample) ~ "F_C2", grepl("^F_C3_", sample) ~ "F_C3",
    grepl("^S_C1_", sample) ~ "S_C1", grepl("^S_C2_", sample) ~ "S_C2", grepl("^S_C3_", sample) ~ "S_C3",
    grepl("^F_T1_", sample) ~ "F_T1", grepl("^F_T2_", sample) ~ "F_T2", grepl("^F_T3_", sample) ~ "F_T3",
    grepl("^S_T1_", sample) ~ "S_T1", grepl("^S_T2_", sample) ~ "S_T2", grepl("^S_T3_", sample) ~ "S_T3",
  )) %>% group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "group")
sum(cc_ont_fl18s_genus_assigned_avg)

cc_ont_fl18s_genus_assigned_avg = cc_ont_fl18s_genus_assigned_avg %>% 
  t() %>% data.frame() %>% rownames_to_column(., var = "Group") %>% as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, "F_C1" ~ "CS1-AF", "F_C2" ~ "CS2-AF", "F_C3" ~ "CS3-AF",
                             "S_C1" ~ "CS1-PS", "S_C2" ~ "CS2-PS", "S_C3" ~ "CS3-PS",
                             "F_T1" ~ "ES1-AF", "F_T2" ~ "ES2-AF", "F_T3" ~ "ES3-AF", 
                             "S_T1" ~ "ES1-PS", "S_T2" ~ "ES2-PS", "S_T3" ~ "ES3-PS")) %>% arrange(Sample)
sum(cc_ont_fl18s_genus_assigned_avg[,"Proportion"])

cc_ont_fl18s_genus_assigned_avg$Sample = factor(cc_ont_fl18s_genus_assigned_avg$Sample, levels = c("CS1-AF", "CS2-AF", "CS3-AF",
                                                                                                   "CS1-PS", "CS2-PS", "CS3-PS",
                                                                                                   "ES1-AF", "ES2-AF", "ES3-AF",
                                                                                                   "ES1-PS", "ES2-PS", "ES3-PS")) 

cc_ont_fl18s_genus_assigned_avg$Group = factor(cc_ont_fl18s_genus_assigned_avg$Group, levels = c("Thalassiosira", "MAST-3J", "Skeletonema", "Cyclotella", "Cylindrotheca",          
                                                                                                 "Tripos", "Navicula", "Chaetoceros", "Ephelota", "Ostreococcus",           
                                                                                                 "Prorocentrum", "Pelagostrobilidium", "Teleaulax", "Eucampia", "Spirotontonia",
                                                                                                 "Others"))
str(cc_ont_fl18s_genus_assigned_avg)







# Combine Datasets
# 16S Phylum
cc_illu_16sv4_phylum_assigned_avg$Seq = "ILLUMINA"
cc_ont_16sv4_phylum_assigned_avg$Seq = "ONT_SHORT"
cc_ont_fl16s_phylum_assigned_avg$Seq = "ONT_LONG"

cc_phylum_combined_16s_avg = bind_rows(
  cc_illu_16sv4_phylum_assigned_avg, 
  cc_ont_16sv4_phylum_assigned_avg,
  cc_ont_fl16s_phylum_assigned_avg
) %>%
  mutate(
    Seq = factor(Seq, levels = c("ILLUMINA", "ONT_SHORT", "ONT_LONG"))) %>%
  arrange(Sample, Seq) %>%
  mutate(Plot_X = paste(Sample, Seq, sep = "---"),
         Flow_Group = paste(Sample, Group, sep = "_")
  ) %>%
  mutate(Plot_X = factor(Plot_X, levels = unique(Plot_X)))

cc_phylum_16s_avg = rbind(cc_illu_16sv4_phylum_assigned_avg, cc_ont_16sv4_phylum_assigned_avg, cc_ont_fl16s_phylum_assigned_avg)
cc_phylum_16s_avg$Seq = factor (cc_phylum_16s_avg$Seq, levels = c( "ILLUMINA", "ONT_SHORT", "ONT_LONG"))
cc_phylum_16s_avg = cc_phylum_16s_avg %>%
  pivot_wider(
    names_from = Seq,
    values_from = Proportion,
    values_fill = 0
  )
str(cc_phylum_16s_avg)


# 18S Subdivision
cc_illu_18sv9_subdivision_assigned_avg$Seq = "ILLUMINA"
cc_ont_18sv9_subdivision_assigned_avg$Seq = "ONT_SHORT"
cc_ont_fl18s_subdivision_assigned_avg$Seq = "ONT_LONG"

cc_subdivision_combined_18s_avg = bind_rows(
  cc_illu_18sv9_subdivision_assigned_avg, 
  cc_ont_18sv9_subdivision_assigned_avg,
  cc_ont_fl18s_subdivision_assigned_avg
) %>%
  mutate(
    Seq = factor(Seq, levels = c("ILLUMINA", "ONT_SHORT", "ONT_LONG"))) %>%
  arrange(Sample, Seq) %>%
  mutate(Plot_X = paste(Sample, Seq, sep = "---"),
         Flow_Group = paste(Sample, Group, sep = "_")
  ) %>%
  mutate(Plot_X = factor(Plot_X, levels = unique(Plot_X)))


cc_subdivision_18s_avg = rbind(cc_illu_18sv9_subdivision_assigned_avg, cc_ont_18sv9_subdivision_assigned_avg, cc_ont_fl18s_subdivision_assigned_avg)
cc_subdivision_18s_avg$Seq = factor (cc_subdivision_18s_avg$Seq, levels = c( "ILLUMINA", "ONT_SHORT", "ONT_LONG"))
cc_subdivision_18s_avg = cc_subdivision_18s_avg %>%
  pivot_wider(
    names_from = Seq,
    values_from = Proportion,
    values_fill = 0
  )
str(cc_subdivision_18s_avg)


# 16S Genus
cc_illu_16sv4_genus_assigned_avg$Seq = "ILLUMINA"
cc_ont_16sv4_genus_assigned_avg$Seq = "ONT_SHORT"
cc_ont_fl16s_genus_assigned_avg$Seq = "ONT_LONG"

cc_genus_combined_16s_avg = bind_rows(
  cc_illu_16sv4_genus_assigned_avg, 
  cc_ont_16sv4_genus_assigned_avg,
  cc_ont_fl16s_genus_assigned_avg
) %>%
  mutate(
    Seq = factor(Seq, levels = c("ILLUMINA", "ONT_SHORT", "ONT_LONG"))) %>%
  arrange(Sample, Seq) %>%
  mutate(Plot_X = paste(Sample, Seq, sep = "---"),
         Flow_Group = paste(Sample, Group, sep = "_")
  ) %>%
  mutate(Plot_X = factor(Plot_X, levels = unique(Plot_X)))

cc_genus_16s_avg = rbind(cc_illu_16sv4_genus_assigned_avg, cc_ont_16sv4_genus_assigned_avg, cc_ont_fl16s_genus_assigned_avg)
cc_genus_16s_avg$Seq = factor (cc_genus_16s_avg$Seq, levels = c( "ILLUMINA", "ONT_SHORT", "ONT_LONG"))
cc_genus_16s_avg = cc_genus_16s_avg %>%
  pivot_wider(
    names_from = Seq,
    values_from = Proportion,
    values_fill = 0
  )

# 18S Genus
cc_illu_18sv9_genus_assigned_avg$Seq = "ILLUMINA"
cc_ont_18sv9_genus_assigned_avg$Seq = "ONT_SHORT"
cc_ont_fl18s_genus_assigned_avg$Seq = "ONT_LONG"

cc_genus_combined_18s_avg = bind_rows(
  cc_illu_18sv9_genus_assigned_avg, 
  cc_ont_18sv9_genus_assigned_avg,
  cc_ont_fl18s_genus_assigned_avg
) %>%
  mutate(
    Seq = factor(Seq, levels = c("ILLUMINA", "ONT_SHORT", "ONT_LONG"))) %>%
  arrange(Sample, Seq) %>%
  mutate(Plot_X = paste(Sample, Seq, sep = "---"),
         Flow_Group = paste(Sample, Group, sep = "_")
  ) %>%
  mutate(Plot_X = factor(Plot_X, levels = unique(Plot_X)))

cc_genus_18s_avg = rbind(cc_illu_18sv9_genus_assigned_avg, cc_ont_18sv9_genus_assigned_avg, cc_ont_fl18s_genus_assigned_avg)
cc_genus_18s_avg$Seq = factor (cc_genus_18s_avg$Seq, levels = c( "ILLUMINA", "ONT_SHORT", "ONT_LONG"))
cc_genus_18s_avg = cc_genus_18s_avg %>%
  pivot_wider(
    names_from = Seq,
    values_from = Proportion,
    values_fill = 0
  )



## stats
# 16S phylum 
cor(cc_illu_16sv4_phylum_assigned_avg$Proportion, cc_ont_16sv4_phylum_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_16sv4_phylum_assigned_avg$Proportion, cc_ont_16sv4_phylum_assigned_avg$Proportion, method = "pearson") 

cor(cc_ont_16sv4_phylum_assigned_avg$Proportion, cc_ont_fl16s_phylum_assigned_avg$Proportion, method = "pearson")
cor.test(cc_ont_16sv4_phylum_assigned_avg$Proportion, cc_ont_fl16s_phylum_assigned_avg$Proportion, method = "pearson")

cor(cc_illu_16sv4_phylum_assigned_avg$Proportion, cc_ont_fl16s_phylum_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_16sv4_phylum_assigned_avg$Proportion, cc_ont_fl16s_phylum_assigned_avg$Proportion, method = "pearson") 

# 16S genus
cor(cc_illu_16sv4_genus_assigned_avg$Proportion, cc_ont_16sv4_genus_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_16sv4_genus_assigned_avg$Proportion, cc_ont_16sv4_genus_assigned_avg$Proportion, method = "pearson") 

cor(cc_ont_16sv4_genus_assigned_avg$Proportion, cc_ont_fl16s_genus_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_ont_16sv4_genus_assigned_avg$Proportion, cc_ont_fl16s_genus_assigned_avg$Proportion, method = "pearson") 

cor(cc_illu_16sv4_genus_assigned_avg$Proportion, cc_ont_fl16s_genus_assigned_avg$Proportion, method = "pearson")
cor.test(cc_illu_16sv4_genus_assigned_avg$Proportion, cc_ont_fl16s_genus_assigned_avg$Proportion, method = "pearson") 

# 18S division
cor(cc_illu_18sv9_subdivision_assigned_avg$Proportion, cc_ont_18sv9_subdivision_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_18sv9_subdivision_assigned_avg$Proportion, cc_ont_18sv9_subdivision_assigned_avg$Proportion, method = "pearson") 

cor(cc_ont_18sv9_subdivision_assigned_avg$Proportion, cc_ont_fl18s_subdivision_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_ont_18sv9_subdivision_assigned_avg$Proportion, cc_ont_fl18s_subdivision_assigned_avg$Proportion, method = "pearson") 

cor(cc_illu_18sv9_subdivision_assigned_avg$Proportion, cc_ont_fl18s_subdivision_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_18sv9_subdivision_assigned_avg$Proportion, cc_ont_fl18s_subdivision_assigned_avg$Proportion, method = "pearson") 

# 18S genus
cor(cc_illu_18sv9_genus_assigned_avg$Proportion, cc_ont_18sv9_genus_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_18sv9_genus_assigned_avg$Proportion, cc_ont_18sv9_genus_assigned_avg$Proportion, method = "pearson") 

cor(cc_ont_18sv9_genus_assigned_avg$Proportion, cc_ont_fl18s_genus_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_ont_18sv9_genus_assigned_avg$Proportion, cc_ont_fl18s_genus_assigned_avg$Proportion, method = "pearson") 

cor(cc_illu_18sv9_genus_assigned_avg$Proportion, cc_ont_fl18s_genus_assigned_avg$Proportion, method = "pearson") 
cor.test(cc_illu_18sv9_genus_assigned_avg$Proportion, cc_ont_fl18s_genus_assigned_avg$Proportion, method = "pearson") 


summary(abs(cc_illu_16sv4_phylum_assigned_avg$Proportion - cc_ont_16sv4_phylum_assigned_avg$Proportion)) 
summary(abs(cc_ont_16sv4_phylum_assigned_avg$Proportion - cc_ont_fl16s_phylum_assigned_avg$Proportion)) 
summary(abs(cc_illu_16sv4_phylum_assigned_avg$Proportion - cc_ont_fl16s_phylum_assigned_avg$Proportion)) 

summary(abs(cc_illu_16sv4_genus_assigned_avg$Proportion - cc_ont_16sv4_genus_assigned_avg$Proportion)) 
summary(abs(cc_ont_16sv4_genus_assigned_avg$Proportion - cc_ont_fl16s_genus_assigned_avg$Proportion))  
summary(abs(cc_illu_16sv4_genus_assigned_avg$Proportion - cc_ont_fl16s_genus_assigned_avg$Proportion)) 

summary(abs(cc_illu_18sv9_subdivision_assigned_avg$Proportion - cc_ont_18sv9_subdivision_assigned_avg$Proportion)) 
summary(abs(cc_ont_18sv9_subdivision_assigned_avg$Proportion - cc_ont_fl18s_subdivision_assigned_avg$Proportion)) 
summary(abs(cc_illu_18sv9_subdivision_assigned_avg$Proportion - cc_ont_fl18s_subdivision_assigned_avg$Proportion)) 

summary(abs(cc_illu_18sv9_genus_assigned_avg$Proportion - cc_ont_18sv9_genus_assigned_avg$Proportion)) 
summary(abs(cc_ont_18sv9_genus_assigned_avg$Proportion - cc_ont_fl18s_genus_assigned_avg$Proportion)) 
summary(abs(cc_illu_18sv9_genus_assigned_avg$Proportion - cc_ont_fl18s_genus_assigned_avg$Proportion)) 







## PLOT
# Scatter Plot
# 16s phylum 
cc_plot_16s_phylum_illu_ontshort = ggplot(cc_phylum_16s_avg, aes(x = ILLUMINA, y = ONT_SHORT)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Pseudomonadota" = "#aec7e8", "Bacteroidota" = "#ffbb78", "Cyanobacteriota" = "#C8E6C9", "Actinomycetota" = "#ff9896", "Planctomycetota" = "#c5b0d5",
                                "Verrucomicrobiota" = "#c49c94", "Marinisomatota" = "#f7b6d2", "Myxococcota_A_473307" = "#c7c7c7", "SAR324" = "#dbdb8d", "Bacillota_I" = "#9edae5",
                                "Fusobacteriota" = "#e5c494", "Bdellovibrionota_473306" = "#b3e2cd", "Campylobacterota_A"= "#fff2cc", "Thermosulfidibacterota" = "#bcbddc", "Desulfobacterota_G_459543" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_16s_phylum_ontshort_ontlong = ggplot(cc_phylum_16s_avg, aes(x = ONT_SHORT, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Pseudomonadota" = "#aec7e8", "Bacteroidota" = "#ffbb78", "Cyanobacteriota" = "#C8E6C9", "Actinomycetota" = "#ff9896", "Planctomycetota" = "#c5b0d5",
                                "Verrucomicrobiota" = "#c49c94", "Marinisomatota" = "#f7b6d2", "Myxococcota_A_473307" = "#c7c7c7", "SAR324" = "#dbdb8d", "Bacillota_I" = "#9edae5",
                                "Fusobacteriota" = "#e5c494", "Bdellovibrionota_473306" = "#b3e2cd", "Campylobacterota_A"= "#fff2cc", "Thermosulfidibacterota" = "#bcbddc", "Desulfobacterota_G_459543" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_16s_phylum_illu_ontlong = ggplot(cc_phylum_16s_avg, aes(x = ILLUMINA, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Pseudomonadota" = "#aec7e8", "Bacteroidota" = "#ffbb78", "Cyanobacteriota" = "#C8E6C9", "Actinomycetota" = "#ff9896", "Planctomycetota" = "#c5b0d5",
                                "Verrucomicrobiota" = "#c49c94", "Marinisomatota" = "#f7b6d2", "Myxococcota_A_473307" = "#c7c7c7", "SAR324" = "#dbdb8d", "Bacillota_I" = "#9edae5",
                                "Fusobacteriota" = "#e5c494", "Bdellovibrionota_473306" = "#b3e2cd", "Campylobacterota_A"= "#fff2cc", "Thermosulfidibacterota" = "#bcbddc", "Desulfobacterota_G_459543" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

# 18s subdivision
cc_plot_18s_subdivision_illu_ontshort = ggplot(cc_subdivision_18s_avg, aes(x = ILLUMINA, y = ONT_SHORT)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Gyrista" = "#aec7e8", "Dinoflagellata" = "#ffbb78", "Ciliophora" = "#C8E6C9", "Bigyra" = "#ff9896", "Chlorophyta_X" = "#c5b0d5",
                                "Cryptophyta_X" = "#c49c94", "Cercozoa" = "#f7b6d2", "Haptophyta_X" = "#c7c7c7", "Picozoa_X" = "#dbdb8d", "Rhodophyta_X" = "#9edae5",
                                "Radiolaria" = "#e5c494", "Fungi" = "#b3e2cd", "Ancyromonadida_X"= "#fff2cc", "Apusomonada_X" = "#bcbddc", "Kathablepharida" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_18s_subdivision_ontshort_ontlong = ggplot(cc_subdivision_18s_avg, aes(x = ONT_SHORT, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Gyrista" = "#aec7e8", "Dinoflagellata" = "#ffbb78", "Ciliophora" = "#C8E6C9", "Bigyra" = "#ff9896", "Chlorophyta_X" = "#c5b0d5",
                                "Cryptophyta_X" = "#c49c94", "Cercozoa" = "#f7b6d2", "Haptophyta_X" = "#c7c7c7", "Picozoa_X" = "#dbdb8d", "Rhodophyta_X" = "#9edae5",
                                "Radiolaria" = "#e5c494", "Fungi" = "#b3e2cd", "Ancyromonadida_X"= "#fff2cc", "Apusomonada_X" = "#bcbddc", "Kathablepharida" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_18s_subdivision_illu_ontlong = ggplot(cc_subdivision_18s_avg, aes(x = ILLUMINA, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Gyrista" = "#aec7e8", "Dinoflagellata" = "#ffbb78", "Ciliophora" = "#C8E6C9", "Bigyra" = "#ff9896", "Chlorophyta_X" = "#c5b0d5",
                                "Cryptophyta_X" = "#c49c94", "Cercozoa" = "#f7b6d2", "Haptophyta_X" = "#c7c7c7", "Picozoa_X" = "#dbdb8d", "Rhodophyta_X" = "#9edae5",
                                "Radiolaria" = "#e5c494", "Fungi" = "#b3e2cd", "Ancyromonadida_X"= "#fff2cc", "Apusomonada_X" = "#bcbddc", "Kathablepharida" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))




# 16S genus
cc_plot_16s_genus_illu_ontshort = ggplot(cc_genus_16s_avg, aes(x = ILLUMINA, y = ONT_SHORT)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Marinagarivorans" = "#aec7e8", "Actinomarina" = "#ffbb78", "Parasynechococcus" = "#C8E6C9", "HIMB59" = "#ff9896", "Cellvibrio" = "#c5b0d5",
                                "Luminiphilus" = "#c49c94", "MED-G82" = "#f7b6d2", "Pelagibacter_A_533952" = "#c7c7c7", "TMED14" = "#dbdb8d", "UBA10364" = "#9edae5",
                                "UBA1275" = "#e5c494", "BACL14" = "#b3e2cd", "UBA1014" = "#fff2cc", "UBA3031" = "#bcbddc", "RS62" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_16s_genus_ontshort_ontlong = ggplot(cc_genus_16s_avg, aes(x = ONT_SHORT, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Marinagarivorans" = "#aec7e8", "Actinomarina" = "#ffbb78", "Parasynechococcus" = "#C8E6C9", "HIMB59" = "#ff9896", "Cellvibrio" = "#c5b0d5",
                                "Luminiphilus" = "#c49c94", "MED-G82" = "#f7b6d2", "Pelagibacter_A_533952" = "#c7c7c7", "TMED14" = "#dbdb8d", "UBA10364" = "#9edae5",
                                "UBA1275" = "#e5c494", "BACL14" = "#b3e2cd", "UBA1014" = "#fff2cc", "UBA3031" = "#bcbddc", "RS62" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_16s_genus_illu_ontlong = ggplot(cc_genus_16s_avg, aes(x = ILLUMINA, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Marinagarivorans" = "#aec7e8", "Actinomarina" = "#ffbb78", "Parasynechococcus" = "#C8E6C9", "HIMB59" = "#ff9896", "Cellvibrio" = "#c5b0d5",
                                "Luminiphilus" = "#c49c94", "MED-G82" = "#f7b6d2", "Pelagibacter_A_533952" = "#c7c7c7", "TMED14" = "#dbdb8d", "UBA10364" = "#9edae5",
                                "UBA1275" = "#e5c494", "BACL14" = "#b3e2cd", "UBA1014" = "#fff2cc", "UBA3031" = "#bcbddc", "RS62" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))




# 18S genus
cc_plot_18s_genus_illu_ontshort = ggplot(cc_genus_18s_avg, aes(x = ILLUMINA, y = ONT_SHORT)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Thalassiosira" = "#aec7e8", "MAST-3J" = "#ffbb78", "Skeletonema" = "#C8E6C9", "Cyclotella" = "#ff9896", "Cylindrotheca" = "#c5b0d5",
                                "Tripos" = "#c49c94", "Navicula" = "#f7b6d2", "Chaetoceros" = "#c7c7c7", "Ephelota" = "#dbdb8d", "Ostreococcus" = "#9edae5",
                                "Prorocentrum" = "#e5c494", "Pelagostrobilidium" = "#b3e2cd", "Teleaulax" = "#fff2cc", "Eucampia" = "#bcbddc", "Spirotontonia" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_18s_genus_ontshort_ontlong = ggplot(cc_genus_18s_avg, aes(x = ONT_SHORT, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Thalassiosira" = "#aec7e8", "MAST-3J" = "#ffbb78", "Skeletonema" = "#C8E6C9", "Cyclotella" = "#ff9896", "Cylindrotheca" = "#c5b0d5",
                                "Tripos" = "#c49c94", "Navicula" = "#f7b6d2", "Chaetoceros" = "#c7c7c7", "Ephelota" = "#dbdb8d", "Ostreococcus" = "#9edae5",
                                "Prorocentrum" = "#e5c494", "Pelagostrobilidium" = "#b3e2cd", "Teleaulax" = "#fff2cc", "Eucampia" = "#bcbddc", "Spirotontonia" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_plot_18s_genus_illu_ontlong = ggplot(cc_genus_18s_avg, aes(x = ILLUMINA, y = ONT_LONG)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(color = Group), size = 2, shape = 16, alpha = 1, show.legend = F) +
  scale_color_manual(values = c("Thalassiosira" = "#aec7e8", "MAST-3J" = "#ffbb78", "Skeletonema" = "#C8E6C9", "Cyclotella" = "#ff9896", "Cylindrotheca" = "#c5b0d5",
                                "Tripos" = "#c49c94", "Navicula" = "#f7b6d2", "Chaetoceros" = "#c7c7c7", "Ephelota" = "#dbdb8d", "Ostreococcus" = "#9edae5",
                                "Prorocentrum" = "#e5c494", "Pelagostrobilidium" = "#b3e2cd", "Teleaulax" = "#fff2cc", "Eucampia" = "#bcbddc", "Spirotontonia" = "#fdcdac",
                                "Others" = "#e0e0e080")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("")+ xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x =element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))


# Bar Plot
library(ggalluvial)
cc_panel_combined_16s_phylum_avg = ggplot(cc_phylum_combined_16s_avg, aes(x = Plot_X, y = Proportion, alluvium = Flow_Group, stratum = Group)) +
  geom_alluvium(aes(fill = Group), alpha = 0.5, color = NA, width = 2/3) +
  geom_stratum(aes(fill = Group), alpha = 1, color = NA, width = 2/3) +
  geom_vline(xintercept = seq(3.5, (12 - 1) * 3 + 0.5, by = 3), linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Pseudomonadota" = "#aec7e8", "Bacteroidota" = "#ffbb78", "Cyanobacteriota" = "#C8E6C9", "Actinomycetota" = "#ff9896", "Planctomycetota" = "#c5b0d5",
                               "Verrucomicrobiota" = "#c49c94", "Marinisomatota" = "#f7b6d2", "Myxococcota_A_473307" = "#c7c7c7", "SAR324" = "#dbdb8d", "Bacillota_I" = "#9edae5",
                               "Fusobacteriota" = "#e5c494", "Bdellovibrionota_473306" = "#b3e2cd", "Campylobacterota_A"= "#fff2cc", "Thermosulfidibacterota" = "#bcbddc", "Desulfobacterota_G_459543" = "#fdcdac",
                               "Others" = "#e0e0e080")) +
  scale_y_continuous(limits = c(0, 100.1), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_panel_combined_16s_genus_avg = ggplot(cc_genus_combined_16s_avg, aes(x = Plot_X, y = Proportion, alluvium = Flow_Group, stratum = Group)) +
  geom_alluvium(aes(fill = Group), alpha = 0.5, color = NA, width = 2/3) +
  geom_stratum(aes(fill = Group), alpha = 1, color = NA, width = 2/3) +
  geom_vline(xintercept = seq(3.5, (12 - 1) * 3 + 0.5, by = 3), linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Marinagarivorans" = "#aec7e8", "Actinomarina" = "#ffbb78", "Parasynechococcus" = "#C8E6C9", "HIMB59" = "#ff9896", "Cellvibrio" = "#c5b0d5",
                               "Luminiphilus" = "#c49c94", "MED-G82" = "#f7b6d2", "Pelagibacter_A_533952" = "#c7c7c7", "TMED14" = "#dbdb8d", "UBA10364" = "#9edae5",
                               "UBA1275" = "#e5c494", "BACL14" = "#b3e2cd", "UBA1014" = "#fff2cc", "UBA3031" = "#bcbddc", "RS62" = "#fdcdac",
                               "Others" = "#e0e0e080")) +
  scale_y_continuous(limits = c(0, 100.1), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_panel_combined_18s_subdivision_avg = ggplot(cc_subdivision_combined_18s_avg, aes(x = Plot_X, y = Proportion, alluvium = Flow_Group, stratum = Group)) +
  geom_alluvium(aes(fill = Group), alpha = 0.5, color = NA, width = 2/3) +
  geom_stratum(aes(fill = Group), alpha = 1, color = NA, width = 2/3) +
  geom_vline(xintercept = seq(3.5, (12 - 1) * 3 + 0.5, by = 3), linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Gyrista" = "#aec7e8", "Dinoflagellata" = "#ffbb78", "Ciliophora" = "#C8E6C9", "Bigyra" = "#ff9896", "Chlorophyta_X" = "#c5b0d5",
                               "Cryptophyta_X" = "#c49c94", "Cercozoa" = "#f7b6d2", "Haptophyta_X" = "#c7c7c7", "Picozoa_X" = "#dbdb8d", "Rhodophyta_X" = "#9edae5",
                               "Radiolaria" = "#e5c494", "Fungi" = "#b3e2cd", "Ancyromonadida_X"= "#fff2cc", "Apusomonada_X" = "#bcbddc", "Kathablepharida" = "#fdcdac",
                               "Others" = "#e0e0e080")) +
  scale_y_continuous(limits = c(0, 100.1), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

cc_panel_combined_18s_genus_avg = ggplot(cc_genus_combined_18s_avg, aes(x = Plot_X, y = Proportion, alluvium = Flow_Group, stratum = Group)) +
  geom_alluvium(aes(fill = Group), alpha = 0.5, color = NA, width = 2/3) +
  geom_stratum(aes(fill = Group), alpha = 1, color = NA, width = 2/3) +
  geom_vline(xintercept = seq(3.5, (12 - 1) * 3 + 0.5, by = 3), linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Thalassiosira" = "#aec7e8", "MAST-3J" = "#ffbb78", "Skeletonema" = "#C8E6C9", "Cyclotella" = "#ff9896", "Cylindrotheca" = "#c5b0d5",
                               "Tripos" = "#c49c94", "Navicula" = "#f7b6d2", "Chaetoceros" = "#c7c7c7", "Ephelota" = "#dbdb8d", "Ostreococcus" = "#9edae5",
                               "Prorocentrum" = "#e5c494", "Pelagostrobilidium" = "#b3e2cd", "Teleaulax" = "#fff2cc", "Eucampia" = "#bcbddc", "Spirotontonia" = "#fdcdac",
                               "Others" = "#e0e0e080")) +
  scale_y_continuous(limits = c(0, 100.1), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))


## Figure 
library(patchwork)

figure_cc_scatter_plot = cc_plot_16s_phylum_illu_ontshort + cc_plot_16s_phylum_ontshort_ontlong + cc_plot_16s_phylum_illu_ontlong + 
                         cc_plot_16s_genus_illu_ontshort + cc_plot_16s_genus_ontshort_ontlong + cc_plot_16s_genus_illu_ontlong + 
                         cc_plot_18s_subdivision_illu_ontshort + cc_plot_18s_subdivision_ontshort_ontlong + cc_plot_18s_subdivision_illu_ontlong +
                         cc_plot_18s_genus_illu_ontshort + cc_plot_18s_genus_ontshort_ontlong + cc_plot_18s_genus_illu_ontlong +  
                         plot_layout(ncol = 3, nrow = 4, byrow = T)

ggsave("figure_cc_scatter_plot.pdf", plot = figure_cc_scatter_plot, width = 17, height = 20, units = "cm", device = "pdf")


figure_cc_bar_plot = cc_panel_combined_16s_phylum_avg + cc_panel_combined_16s_genus_avg +
                     cc_panel_combined_18s_subdivision_avg + cc_panel_combined_18s_genus_avg +
                     plot_layout(ncol = 1, nrow = 4, byrow = T)

ggsave("figure_cc_bar_plot.pdf", plot = figure_cc_bar_plot, width = 16, height = 20, units = "cm", device = "pdf")























### NMDS --------------------------------------------------------------------------------
library(tidyverse)
library(vegan)

## data
bray_illu_16sv4 = vegdist(otu_illu_16sv4_rarefied, method = "bray")
bray_ont_16sv4 = vegdist(otu_ont_16sv4_rarefied, method = "bray")
bray_ont_fl16s = vegdist(otu_ont_fl16s_rarefied, method = "bray")
bray_illu_18sv9 = vegdist(otu_illu_18sv9_rarefied, method = "bray")
bray_ont_18sv9 = vegdist(otu_ont_18sv9_rarefied, method = "bray")
bray_ont_fl18s = vegdist(otu_ont_fl18s_rarefied, method = "bray")

nmds_illu_16sv4_bray_table = data.frame(NMDS1 = metaMDS(bray_illu_16sv4, k = 2)$points[,1], NMDS2 = metaMDS(bray_illu_16sv4, k = 2) $points[,2])
nmds_illu_16sv4_bray_table = nmds_illu_16sv4_bray_table %>% 
                             mutate(sample = rownames(nmds_illu_16sv4_bray_table)) %>% 
                             mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
                             mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_illu_16sv4_bray_table$environment = as.factor(nmds_illu_16sv4_bray_table$environment)
nmds_illu_16sv4_bray_table$group = as.factor(nmds_illu_16sv4_bray_table$group)
str(nmds_illu_16sv4_bray_table)

nmds_ont_16sv4_bray_table = data.frame(NMDS1 = metaMDS(bray_ont_16sv4, k = 2)$points[,1], NMDS2 = metaMDS(bray_ont_16sv4, k = 2) $points[,2])
nmds_ont_16sv4_bray_table = nmds_ont_16sv4_bray_table %>% 
  mutate(sample = rownames(nmds_ont_16sv4_bray_table)) %>% 
  mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_ont_16sv4_bray_table$environment = as.factor(nmds_ont_16sv4_bray_table$environment)
nmds_ont_16sv4_bray_table$group = as.factor(nmds_ont_16sv4_bray_table$group)
str(nmds_ont_16sv4_bray_table)

nmds_ont_fl16s_bray_table = data.frame(NMDS1 = metaMDS(bray_ont_fl16s, k = 2)$points[,1], NMDS2 = metaMDS(bray_ont_fl16s, k = 2) $points[,2])
nmds_ont_fl16s_bray_table = nmds_ont_fl16s_bray_table %>% 
  mutate(sample = rownames(nmds_ont_fl16s_bray_table)) %>% 
  mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_ont_fl16s_bray_table$environment = as.factor(nmds_ont_fl16s_bray_table$environment)
nmds_ont_fl16s_bray_table$group = as.factor(nmds_ont_fl16s_bray_table$group)
str(nmds_ont_fl16s_bray_table)


nmds_illu_18sv9_bray_table = data.frame(NMDS1 = metaMDS(bray_illu_18sv9, k = 2)$points[,1], NMDS2 = metaMDS(bray_illu_18sv9, k = 2) $points[,2])
nmds_illu_18sv9_bray_table = nmds_illu_18sv9_bray_table %>% 
  mutate(sample = rownames(nmds_illu_18sv9_bray_table)) %>% 
  mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_illu_18sv9_bray_table$environment = as.factor(nmds_illu_18sv9_bray_table$environment)
nmds_illu_18sv9_bray_table$group = as.factor(nmds_illu_18sv9_bray_table$group)
str(nmds_illu_18sv9_bray_table)

nmds_ont_18sv9_bray_table = data.frame(NMDS1 = metaMDS(bray_ont_18sv9, k = 2)$points[,1], NMDS2 = metaMDS(bray_ont_18sv9, k = 2) $points[,2])
nmds_ont_18sv9_bray_table = nmds_ont_18sv9_bray_table %>% 
  mutate(sample = rownames(nmds_ont_18sv9_bray_table)) %>% 
  mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_ont_18sv9_bray_table$environment = as.factor(nmds_ont_18sv9_bray_table$environment)
nmds_ont_18sv9_bray_table$group = as.factor(nmds_ont_18sv9_bray_table$group)
str(nmds_ont_18sv9_bray_table)

nmds_ont_fl18s_bray_table = data.frame(NMDS1 = metaMDS(bray_ont_fl18s, k = 2)$points[,1], NMDS2 = metaMDS(bray_ont_fl18s, k = 2) $points[,2])
nmds_ont_fl18s_bray_table = nmds_ont_fl18s_bray_table %>% 
  mutate(sample = rownames(nmds_ont_fl18s_bray_table)) %>% 
  mutate(environment = case_when(grepl("C", sample) ~ "C", grepl("T", sample) ~ "T")) %>% 
  mutate(group = case_when(grepl("^F_C", sample) ~ "F_C", grepl("^F_T", sample) ~ "F_T", grepl("^S_C", sample) ~ "S_C", grepl("^S_T", sample) ~ "S_T"))
nmds_ont_fl18s_bray_table$environment = as.factor(nmds_ont_fl18s_bray_table$environment)
nmds_ont_fl18s_bray_table$group = as.factor(nmds_ont_fl18s_bray_table$group)
str(nmds_ont_fl18s_bray_table)

bray_16s_between_env = c(as.matrix(vegdist(otu_illu_16sv4_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_illu_16sv4_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_illu_16sv4_rarefied, method = "bray"))))],
                         as.matrix(vegdist(otu_ont_16sv4_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_ont_16sv4_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_ont_16sv4_rarefied, method = "bray"))))],
                         as.matrix(vegdist(otu_ont_fl16s_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_ont_fl16s_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_ont_fl16s_rarefied, method = "bray"))))])

bray_16s_between_env = data.frame(distance = bray_16s_between_env, seq = c(rep("ILLUMINA", 324), rep("ONTSHORT", 324), rep("ONTLONG", 324)))
bray_16s_between_env$seq = factor (bray_16s_between_env$seq, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
bray_16s_between_env$dissimilarity = bray_16s_between_env$distance *100
str(bray_16s_between_env)

bray_18s_between_env = c(as.matrix(vegdist(otu_illu_18sv9_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_illu_18sv9_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_illu_18sv9_rarefied, method = "bray"))))],
                         as.matrix(vegdist(otu_ont_18sv9_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_ont_18sv9_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_ont_18sv9_rarefied, method = "bray"))))],
                         as.matrix(vegdist(otu_ont_fl18s_rarefied, method = "bray"))[grepl("C", rownames(as.matrix(vegdist(otu_ont_fl18s_rarefied, method = "bray")))), grepl("T", colnames(as.matrix(vegdist(otu_ont_fl18s_rarefied, method = "bray"))))])
bray_18s_between_env = data.frame(distance = bray_18s_between_env, seq = c(rep("ILLUMINA", 324), rep("ONTSHORT", 324), rep("ONTLONG", 324)))
bray_18s_between_env$dissimilarity = bray_18s_between_env$distance *100
bray_18s_between_env$seq = factor (bray_18s_between_env$seq, levels = c("ILLUMINA","ONTSHORT", "ONTLONG"))
str(bray_18s_between_env)

## stats
metaMDS(bray_illu_16sv4, k = 2)
metaMDS(bray_ont_16sv4, k = 2)
monoMDS(bray_ont_fl16s, k = 2)
metaMDS(bray_illu_18sv9, k = 2)
metaMDS(bray_ont_18sv9, k = 2)
metaMDS(bray_ont_fl18s, k = 2)

betadisper(bray_illu_16sv4, nmds_illu_16sv4_bray_table$environment, type = "median") 
betadisper(bray_ont_16sv4, nmds_ont_16sv4_bray_table$environment, type = "median") 
betadisper(bray_ont_fl16s, nmds_ont_fl16s_bray_table$environment, type = "median")
betadisper(bray_illu_18sv9, nmds_illu_18sv9_bray_table$environment, type = "median") 
betadisper(bray_ont_18sv9, nmds_ont_18sv9_bray_table$environment, type = "median") 
betadisper(bray_ont_fl18s, nmds_ont_fl18s_bray_table$environment, type = "median") 

adonis2(bray_illu_16sv4 ~ nmds_illu_16sv4_bray_table$environment, permutations = 999) 
adonis2(bray_ont_16sv4 ~ nmds_ont_16sv4_bray_table$environment, permutations = 999) 
adonis2(bray_ont_fl16s ~ nmds_ont_fl16s_bray_table$environment, permutations = 999) 
adonis2(bray_illu_18sv9 ~ nmds_illu_18sv9_bray_table$environment, permutations = 999)
adonis2(bray_ont_18sv9 ~ nmds_ont_18sv9_bray_table$environment, permutations = 999) 
adonis2(bray_ont_fl18s ~ nmds_ont_fl18s_bray_table$environment, permutations = 999) 

protest(metaMDS(bray_illu_16sv4, k = 2), metaMDS(bray_ont_16sv4, k = 2), permutations = 999) 
protest(metaMDS(bray_ont_16sv4, k = 2), metaMDS(bray_ont_fl16s, k = 2), permutations = 999) 
protest(metaMDS(bray_illu_16sv4, k = 2), metaMDS(bray_ont_fl16s, k = 2), permutations = 999) 
protest(metaMDS(bray_illu_18sv9, k = 2), metaMDS(bray_ont_18sv9, k = 2), permutations = 999) 
protest(metaMDS(bray_ont_18sv9, k = 2), metaMDS(bray_ont_fl18s, k = 2), permutations = 999) 
protest(metaMDS(bray_illu_18sv9, k = 2), metaMDS(bray_ont_fl18s, k = 2), permutations = 999) 

summary(bray_16s_between_env[bray_16s_between_env$seq == "ILLUMINA", ]$dissimilarity)
sd(bray_16s_between_env[bray_16s_between_env$seq == "ILLUMINA", ]$dissimilarity)
summary(bray_16s_between_env[bray_16s_between_env$seq == "ONTSHORT", ]$dissimilarity)
sd(bray_16s_between_env[bray_16s_between_env$seq == "ONTSHORT", ]$dissimilarity)
summary(bray_16s_between_env[bray_16s_between_env$seq == "ONTLONG", ]$dissimilarity)
sd(bray_16s_between_env[bray_16s_between_env$seq == "ONTLONG", ]$dissimilarity)
summary(bray_18s_between_env[bray_18s_between_env$seq == "ILLUMINA", ]$dissimilarity)
sd(bray_18s_between_env[bray_18s_between_env$seq == "ILLUMINA", ]$dissimilarity)
summary(bray_18s_between_env[bray_18s_between_env$seq == "ONTSHORT", ]$dissimilarity)
sd(bray_18s_between_env[bray_18s_between_env$seq == "ONTSHORT", ]$dissimilarity)
summary(bray_18s_between_env[bray_18s_between_env$seq == "ONTLONG", ]$dissimilarity)
sd(bray_18s_between_env[bray_18s_between_env$seq == "ONTLONG", ]$dissimilarity)


## plot
nmds_figure_illu_16sv4 = ggplot(nmds_illu_16sv4_bray_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-4, 6), breaks = c(-4, 1, 6)) +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

nmds_figure_ont_16sv4 = ggplot(nmds_ont_16sv4_bray_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-6, 8), breaks = c(-6, 1, 8)) +
  scale_y_continuous(limits = c(-5, 5), breaks = c(-5, 0, 5)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

nmds_figure_ont_fl16s = ggplot(nmds_ont_fl16s_bray_table, aes(x = NMDS1, y = NMDS2, shape = environment)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-7, 9), breaks = c(-7, 1, 9)) +
  scale_y_continuous(limits = c(-6, 6), breaks = c(-6, 0, 6)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y =element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 


nmds_figure_illu_18sv9 = ggplot(nmds_illu_18sv9_bray_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

nmds_figure_ont_18sv9 = ggplot(nmds_ont_18sv9_bray_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

nmds_figure_ont_fl18s = ggplot(nmds_ont_fl18s_bray_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = environment), shape = 16, size = 3, alpha = 1) +
  stat_ellipse(aes(color = environment, fill = environment), geom = "polygon", linetype = "solid", linewidth = NA, alpha = 0.4, level = 0.95) +
  scale_color_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_fill_manual(values = c("C" = "#05B9E2", "T" = "#BB9727")) +
  scale_x_continuous(limits = c(-4, 4), breaks = c(-4, 0, 4)) +
  scale_y_continuous(limits = c(-4, 4), breaks = c(-4, 0, 4)) +
  ggtitle("")+ xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size= 8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 



dissimilarity_plot_16s = ggplot(bray_16s_between_env, aes(x = seq, y = dissimilarity)) +
  geom_point(size = 1, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.2, height = 0), aes(color = factor(seq))) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, size = 0.6, aes(color = factor(seq)), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "solid", size = 0.5, color = "black") +
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT\nShort", "ONT\nLong")) +
  scale_y_continuous(limits = c(60, 100), breaks = c(60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("16s") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

dissimilarity_plot_18s = ggplot(bray_18s_between_env, aes(x = seq, y = dissimilarity)) +
  geom_point(size = 1, shape = 16, show.legend = F, alpha = 0.5, position=position_jitter(width = 0.2, height = 0), aes(color = factor(seq))) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, size = 0.6, aes(color = factor(seq)), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "solid", size = 0.5, color = "black") +
  scale_color_manual(values = c("ILLUMINA" = "#A9A9A9", "ONTSHORT" = "#9BC4E4", "ONTLONG" = "#3579A1")) +
  scale_x_discrete(label = c("Illumina","ONT\nShort", "ONT\nLong")) +
  scale_y_continuous(limits = c(80, 100), breaks = c(80, 90, 100)) +
  ggtitle("") + xlab("") + ylab("18s") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))


## Figure 
library(patchwork)
figure_nmds = nmds_figure_illu_16sv4 + nmds_figure_ont_16sv4 + nmds_figure_ont_fl16s + nmds_figure_illu_18sv9 + nmds_figure_ont_18sv9 + nmds_figure_ont_fl18s + plot_layout(ncol = 3, nrow = 2, byrow = T)

ggsave("figure_nmds.pdf", plot = figure_nmds, width = 17, height = 12, units = "cm", device = "pdf")

figure_dissimilarity = dissimilarity_plot_16s + dissimilarity_plot_18s + plot_layout(ncol = 1, nrow = 2, byrow = T)

ggsave("figure_dissimilarity.pdf", plot = figure_dissimilarity, width = 8, height = 10, units = "cm", device = "pdf")



















### dbRDA & HP -------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(vegan)
library(rdacca.hp)
library(patchwork)

## data
rownames(env_factor) = env_factor[,1]
env_factor = env_factor[,-1]
env_factor = env_factor[rep(1:nrow(env_factor), each = 3), ]
env_factor = rbind(env_factor, env_factor)
rownames(env_factor) = rownames(pcoa_illu_16sv4_bray_table)



## stats
dbrda(bray_illu_16sv4 ~ ., env_factor, distance = 'bray') 

dbrda(bray_ont_16sv4 ~ ., env_factor, distance = 'bray')

dbrda(bray_ont_fl16s ~ ., env_factor, distance = 'bray')

dbrda(bray_illu_18sv9 ~ ., env_factor, distance = 'bray')

dbrda(bray_ont_18sv9 ~ ., env_factor, distance = 'bray')

dbrda(bray_ont_fl18s ~ ., env_factor, distance = 'bray')


rdacca.hp(bray_illu_16sv4, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_illu_16sv4, env_factor, method ="dbRDA", type="R2", permutations = 1000)

rdacca.hp(bray_ont_16sv4, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_ont_16sv4, env_factor, method ="dbRDA", type="R2", permutations = 1000)

rdacca.hp(bray_ont_fl16s, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_ont_fl16s, env_factor, method ="dbRDA", type="R2", permutations = 1000)

rdacca.hp(bray_illu_18sv9, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_illu_18sv9, env_factor, method ="dbRDA", type="R2", permutations = 1000)

rdacca.hp(bray_ont_18sv9, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_ont_18sv9, env_factor, method ="dbRDA", type="R2", permutations = 1000)

rdacca.hp(bray_ont_fl18s, env_factor, method ="dbRDA", type="R2", n.perm = 1000)
permu.hp(bray_ont_fl18s, env_factor, method ="dbRDA", type="R2", permutations = 1000)



## plot
hp_plot_illu_16sv4 = ggplot(hp_illu_16sv4, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#A9A9A9", fill = "#A9A9A9") +
  scale_y_continuous(limits = c(0, 14), breaks = c(0, 7, 14)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

hp_plot_ont_16sv4 = ggplot(hp_ont_16sv4, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#9BC4E4", fill = "#9BC4E4") +
  scale_y_continuous(limits = c(0, 14), breaks = c(0, 7, 14)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

hp_plot_ont_fl16s = ggplot(hp_ont_fl16s, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#3579A1", fill = "#3579A1") +
  scale_y_continuous(limits = c(0, 14), breaks = c(0, 7, 14)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 


hp_plot_illu_18sv9 = ggplot(hp_illu_18sv9, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#A9A9A9", fill = "#A9A9A9") +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 9, 18)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

hp_plot_ont_18sv9 = ggplot(hp_ont_18sv9, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#9BC4E4", fill = "#9BC4E4") +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 9, 18)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

hp_plot_ont_fl18s = ggplot(hp_ont_fl18s, aes(x = Factor, y = IndPer)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "#3579A1", fill = "#3579A1") +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 9, 18)) +
  geom_text(aes(label = round(IndPer, 2)), vjust = -0.5, size = 2, family = "Helvetica", fontface = "plain") +
  ggtitle("") + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth = 1, linetype = "solid")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(color = "black", size = 8,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y = element_text(color = "black",size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("hp_plot_ont_fl18s.pdf", plot = hp_plot_ont_fl18s, width = 5, height = 3, units = "cm", device = "pdf")


## Figure
figure_hp = hp_plot_illu_16sv4 + hp_plot_ont_16sv4 + hp_plot_ont_fl16s + hp_plot_illu_18sv9 + hp_plot_ont_18sv9 + hp_plot_ont_fl18s + plot_layout(ncol = 3, nrow = 2, byrow = T)

ggsave("figure_hp.pdf", plot = figure_hp, width = 17.2, height = 9, units = "cm", device = "pdf")





































### MANTEL ---------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(geosphere)
library(linkET)
library(patchwork)

## data
rownames(geo) = geo[,1]
geo = geo[,-1]
geo = geo[rep(1:nrow(geo), each = 3), ]
geo = rbind(geo, geo)
rownames(geo) = rownames(pcoa_illu_16sv4_bray_table)


## stats
mantel.partial(bray_illu_16sv4, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_16sv4, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_16sv4, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_16sv4, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_16sv4, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_16sv4, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 

mantel.partial(bray_ont_16sv4, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_16sv4, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_16sv4, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_16sv4, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_16sv4, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_16sv4, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 

mantel.partial(bray_ont_fl16s, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl16s, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl16s, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl16s, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl16s, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl16s, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 

mantel.partial(bray_illu_18sv9, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_18sv9, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_18sv9, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_18sv9, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_18sv9, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_illu_18sv9, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 

mantel.partial(bray_ont_18sv9, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_18sv9, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_18sv9, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_18sv9, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_18sv9, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_18sv9, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 

mantel.partial(bray_ont_fl18s, dist(env_factor[, "Temperature"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl18s, dist(env_factor[, "Salinity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl18s, dist(env_factor[, "Turbidity"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl18s, dist(env_factor[, "pH"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl18s, dist(env_factor[, "DO"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 
mantel.partial(bray_ont_fl18s, dist(env_factor[, "Chlorophyll"], method = "euclidean"), geosphere::distm(geo[, c('Longitude', 'Laititude')]), method = "pearson", permutations = 999) 


mantel_table_16s = mantel_table_16s %>% 
        mutate(rd = cut(mantel_r, breaks = c(-Inf, 0.25, Inf), labels = c("< 0.25", "≥ 0.25")),
               pd = cut(mantel_p, breaks = c(-Inf, 0.001, Inf), labels = c("< 0.001","≥ 0.001")))

mantel_table_18s = mantel_table_18s %>% 
  mutate(rd = cut(mantel_r, breaks = c(-Inf, 0.25, Inf), labels = c("< 0.25", "≥ 0.25")),
         pd = cut(mantel_p, breaks = c(-Inf, 0.001, Inf), labels = c("< 0.001","≥ 0.001")))


## plot
library(linkET)

mantel_figure_16s = qcorrplot(correlate(env_factor, method = "pearson"), type = "lower", diag = T) +
  geom_square() + 
  geom_couple(aes(colour = pd, size = rd), data = mantel_table_16s, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#05B9E2", "#4FC0EC", "#A2A2A288", "#E6B366", "#BB9727"), limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_colour_manual(values = c("#F27970", "#A2A2A288")) +
  scale_size_manual(values = c(0.5, 1.5)) +
  guides(colour = guide_legend(title = "Mantel's p", order = 1, override.aes = list(size = 3)), size = guide_legend(title = "Mantel's r", order = 2), fill = guide_colorbar(title = "Pearson's r", order = 3, ticks = FALSE)) +
  labs(title = "") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.y.right = element_blank()) +
  theme(plot.title = element_text())

ggsave("mantel_figure_16s.pdf", plot = mantel_figure_16s, width = 10, height = 10, units = "cm", device = "pdf")

mantel_figure_18s = qcorrplot(correlate(env_factor, method = "pearson"), type = "lower", diag = T) +
  geom_square() + 
  geom_couple(aes(colour = pd, size = rd), data = mantel_table_18s, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#05B9E2", "#4FC0EC", "#A2A2A288", "#E6B366", "#BB9727"), limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_colour_manual(values = c("#F27970", "#A2A2A288")) +
  scale_size_manual(values = c(0.5, 1.5)) +
  guides(colour = guide_legend(title = "Mantel's p", order = 1, override.aes = list(size = 3)), size = guide_legend(title = "Mantel's r", order = 2), fill = guide_colorbar(title = "Pearson's r", order = 3, ticks = FALSE)) +
  labs(title = "") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.y.right = element_blank()) +
  theme(plot.title = element_text())

ggsave("mantel_figure_18s.pdf", plot = mantel_figure_18s, width = 10, height = 10, units = "cm", device = "pdf")


## Figure
figure_mantel = mantel_figure_16s + mantel_figure_18s + plot_layout(ncol = 1, nrow = 2, byrow = T)

ggsave("figure_mantel.pdf", plot = figure_mantel, width = 15, height = 15, units = "cm", device = "pdf")