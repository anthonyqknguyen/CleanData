library(edgeR)
library(dplyr)
library(MASS)

ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t", header=TRUE)[,-c(13)]

ko_noz_filt <- ko_noz[rowSums(ko_noz) > 12,]

ko_noz_filt <- filter(ko_noz_filt , rowSums(across(everything(), ~.x==0))<=9)

baseline <- DGEList(counts=ko_noz_filt[,1:6])

treatment <- DGEList(counts=ko_noz_filt[,7:12])


baselineDisp = estimateDisp(baseline)

treatmentDisp = estimateDisp(treatment)

baselineDisp$tagwise.dispersion

baselineParams <- data.frame("gene" = rownames(ko_noz_filt), "disp" = baselineDisp$trended.dispersion, "logCPM" = baselineDisp$AveLogCPM)

treatmentParams <- data.frame("gene" = rownames(ko_noz_filt), "disp" = treatmentDisp$trended.dispersion, "logCPM" = treatmentDisp$AveLogCPM)

theta = (exp(baselineParams$logCPM) + exp(baselineParams$logCPM)^2) / baselineParams$disp

max(rnegbin(exp(baselineParams$logCPM), theta = baselineParams$disp))

max(rnegbin(exp(baselineParams$logCPM), theta = theta))

