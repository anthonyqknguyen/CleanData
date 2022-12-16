
source("Code/RNASeqPower.R")

Data_wd <- "./Data/"


ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t", header=TRUE)[,-c(13)]

ko_noz_filt <- ko_noz[rowSums(ko_noz) > 12,]

ko_noz_filt <- filter(ko_noz_filt , rowSums(across(everything(), ~.x==0))<=9)


snf2_metadata <- read.table(paste(Data_wd, "metadata (1).txt", sep=""), sep="\t", fill = TRUE)

panelist = snf2_metadata$panelist

treatment = snf2_metadata$toothpaste

treatment


results = ofPowerAnalysis(ko_noz_filt, treatment, 2, 5, 10, 1)

plot(rownames(results),rowMeans(results, na.rm=T), main="One Factor Power Analysis on Snf2 Data", xlab = "Number of Replicates", ylab = "Power")

results








