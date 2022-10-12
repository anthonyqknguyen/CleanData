# Setting woking directory
# assumes starting location in terminal is /CleanData
setwd(getwd())


#install.packages("tidyr")
#install.packages("rlang")
library(tidyr)
sessionInfo()
library(rlang)
# Setting the data directory
Data_wd <- "./Data/"

# Reading in the counts table
counts <- read.table(paste(Data_wd, "counts.txt", sep=""), sep="\t", header =TRUE,
                     row.names = 1)


# Reading in the metadata
metadata <- read.delim(paste(Data_wd, "metadata.txt", sep =""), header = TRUE, sep = "\t", row.names = 1)
metadata_filt <- metadata[complete.cases(metadata[ , "Smoker"]),]

#counts_filt <- t(counts)[rownames(t(counts)) %in% rownames(metadata_filt)]

counts <- as.data.frame(t(counts))
#counts_filt <- counts[rownames(counts) %in% rownames(metadata_filt)]

counts_filt <- dplyr::select(counts, rownames(metadata_filt))
# Reading in the taxa data
taxa <- read.table(paste(Data_wd, "taxa.txt", sep=""), sep="\t")

# Reading in Khier data
khier <- read.table(paste(Data_wd, "khier_2020_fixed.txt", sep=""), sep="\t")

# Reading in ko_noz data
ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t")

# Reading in snf2 metadata
snf2_metadata <- read.table(paste(Data_wd, "metadata (1).txt", sep=""), sep="\t", fill = TRUE)

#View(taxa)
#View(metadata)
library("DESeq2")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = metadata_filt, 
                               design = ~Drinks.Coffee)
DESeq(dds, fitType = 'glmGamPoi')
#View(dds$Smoker)
