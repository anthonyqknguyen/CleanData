# Setting woking directory
# assumes starting location in terminal is /CleanData
setwd(getwd())

# Setting the data directory
Data_wd <- "./Data/"

# Reading in the counts table
counts <- read.table(paste(Data_wd, "counts.txt", sep=""), sep="\t")

# Reading in the citizen metadata
citizen_metadata <- read.table(paste(Data_wd, "metadata.txt", sep=""), sep="\t", fill = TRUE)

# Reading in the taxa data
taxa <- read.table(paste(Data_wd, "taxa.txt", sep=""), sep="\t")

# Reading in Khier data
khier <- read.table(paste(Data_wd, "khier_2020_fixed.txt", sep=""), sep="\t")

# Reading in ko_noz data
ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t")

# Reading in snf2 metadata
snf2_metadata <- read.table(paste(Data_wd, "metadata (1).txt", sep=""), sep="\t", fill = TRUE)