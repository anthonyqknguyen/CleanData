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

simData = rnegbin(exp(baselineParams$logCPM), theta = baselineParams$disp)

max(rnegbin(exp(baselineParams$logCPM), theta = theta))

boxplot(baselineDisp$trended.dispersion, treatmentDisp$trended.dispersion, 
        ylab="Distribution of Fitted Dispersion", xlab="Study (Before/After Treatment)", main="Fitted Dispersion Values Before and After Treatment", col=c("firebrick4", "dodgerblue4"))
legend("topright", legend = c("Before Treatment", "After Treatment"), fill = c("firebrick4", "dodgerblue4"))

boxplot(baselineDisp$AveLogCPM, treatmentDisp$AveLogCPM, 
        ylab="Distribution of Average Log CPM (counts per million)", xlab="Study (Before/After Treatment)", main="Fitted Average Log CPM Values Before and After Treatment", col=c("firebrick4", "dodgerblue4"))
legend("topright", legend = c("Before Treatment", "After Treatment"), fill = c("firebrick4", "dodgerblue4"))

plot(sort(log(simData)))

plot(sort(simData), ylim = c(0,100000))





