# Setting woking directory
# assumes starting location in terminal is /CleanData
setwd("C:/Users/antho/OneDrive - Virginia Tech/School/Year 4/Fall_2022/CMDA_4864/CleanData")

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


library(dplyr)

make_hist <- function(data, title= "untitled", save= FALSE, file_name= "untitled.svg") {
  hist.data <- hist(as.numeric(do.call('c', data)), plot=FALSE)
  hist.data$counts = log10(hist.data$counts + 1)
  if (save) {
    svg(file_name)
  }
  plot(hist.data,
       xlab="Counts", ylab="Frequency (log10)", main=title, 
       col="lightblue")
  if (save) {
    dev.off()
  }
}

counts_num <- select(counts, -V1)
ko_num <- select(ko_noz, -ko)
ko_num_before <- ko_num[,1:6]
ko_num_after <- ko_num[,7:12]

make_hist(counts_num, "Citizens Data")
make_hist(ko_num_before, "SnF2 Before Treatment Data")
make_hist(ko_num_after, "SnF2 After Treatment Data")
