# Setting woking directory
# assumes starting location in terminal is /CleanData
setwd(getwd())

# Setting the data directory
Data_wd <- "../Data/"

# Reading in the counts table
counts <- read.table(paste(Data_wd, "counts.txt", sep=""), sep="\t")

# Reading in the metadata
metadata <- read.table(paste(Data_wd, "metadata.txt", sep=""), sep="\t", fill = TRUE)

# Reading in the taxa data
taxa <- read.table(paste(Data_wd, "taxa.txt", sep=""), sep="\t")