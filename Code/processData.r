setwd(getwd())
Data_wd <- "../Data/"
counts <- read.table(paste(Data_wd, "counts.txt", sep=""), sep="\t")
metadata <- read.table(paste(Data_wd, "metadata copy.txt", sep=""), sep="\t")
taxa <- read.table(paste(Data_wd, "taxa.txt", sep=""), sep="\t")
