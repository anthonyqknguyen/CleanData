library(MASS)
library(edgeR)
library(dplyr)

estimatePairedParams = function(counts, panelist, treatment) {
  
  # sort data by treatment
  
  sort_order <- order(treatment)
  
  treatment <- treatment[sort_order]
  
  counts <- counts[, sort_order]
  
  # create the DGEList object using the sorted counts
  
  y = DGEList(counts=counts)
  
  # create design matrix
  
  design = model.matrix(~panelist + treatment)
  
  rownames(design) <- colnames(y)
  
  # estimate dispersions
  
  dispsCR = dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, 
                                          trend = FALSE, AveLogCPM = NULL, min.row.sum = 6, prior.df = 0, span = 0.3, 
                                          grid.npts = 15, grid.range = c(-8, 8))
  
  sample_data <- data.frame(panelist, treatment)
  
  sample_data$libsize <- log(colSums(y$counts))
  
  libsize <- sample_data$libsize
  
  nofit <- c()
  
  fc <- matrix(nrow = dim(y$counts)[1], ncol = length(levels(factor(panelist))) + length(levels(factor(treatment))) - 1)
  
  for (i in 1:dim(y$counts)[1]) {
    f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
    tryCatch({ glm(y$counts[i,] ~ panelist + treatment + 0, offset = libsize, family = f) -> fit },
             warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
    fc[i,] <- fit$coefficients
  }
  
  y <- DGEList(counts = ko_noz_filt[-nofit,])
  
  params = list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
  
  
  # find de for given count data
  
  edgeR_cds <- DGEList(counts = counts[-nofit,])
  
  edgeR_cds <- calcNormFactors(edgeR_cds)
  
  edgeR_cds <- estimateGLMCommonDisp(edgeR_cds, design)
  
  edgeR_cds <- estimateGLMTrendedDisp(edgeR_cds, design)
  
  edgeR_cds <- estimateGLMTagwiseDisp(edgeR_cds, design)
  
  fit <- glmFit(edgeR_cds, design)
  
  res <- glmLRT(fit)$table
  
  pval <- res$PValue
  
  de = pval < 0.05
  
  params[["de"]] <- de
  
  return(params)
  
}


simulatePaired = function(params, n) {
  
  fcCopy = params$fc
  
  fcPanelist = params$fc[, 1:(dim(params$fc)[2] - 1)]
  
  de = params$de
  
  libsize = params$sample_data$libsize
  
  fcCopy[!de] <- 0
  
  
  sim_libsizes <- runif(n * 2, min = min(libsize), max = max(libsize))

  
  panelist_min <- apply(fcPanelist, 1, min)
  
  panelist_max <- apply(fcPanelist, 1, max)
  
  sim_pfc <- t(sapply(1:length(panelist_min), function(i) {
    runif(n, min = panelist_min[i], max = panelist_max[i])
  }))
  
  
  
  m <- matrix(nrow = length(params$dispsCR), ncol = n * 2)
  
  
  fcCopy = fcCopy[, dim(fcCopy)[2]]
  
  for (i in 1:length(params$dispsCR)) {
    for (j in 1:(n * 2)) {
      m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] +
                                       sim_pfc[i, ceiling(j / 2)] +
                                       ifelse(j %% 2 == 0, fcCopy[i], 0)), size = 1 / params$dispsCR[i])
    }
  }
  
  
  
  Panelist <- sort(c(1:n, 1:n))
  
  Treatment <- rep(c("N", "T"), n)
  
  colnames(m) <- paste0("s", Panelist, Treatment)
  
  rownames(m) <- paste0("g", 1:length(params$dispsCR))
  
  return(m)
  
}





evalSimulatedData = function(m, n) {
  
  pval_list_sim = list()
    
  panelistFactor <- factor(sort(c(1:n, 1:n)))
    
  treatmentFactor <- factor(rep(c("N", "T"), n))
    
  design <- model.matrix(~panelistFactor + treatmentFactor)
    
  edgeR_cds <- DGEList(m)
    
  edgeR_cds <- calcNormFactors(edgeR_cds)
    
  edgeR_cds <- estimateGLMCommonDisp(edgeR_cds, design)
    
  edgeR_cds <- estimateGLMTrendedDisp(edgeR_cds, design)
    
  edgeR_cds <- estimateGLMTagwiseDisp(edgeR_cds, design)
    
  fit <- glmFit(edgeR_cds, design)
    
  res <- glmLRT(fit)$table
    
  pvalSim <- res$PValue
    
  padjSim <- p.adjust(pvalSim, method = "BH")
    
  resSim <- cbind(pvalSim, padjSim)
    
  pval_list_sim[["er"]] <- as.matrix(resSim)
  
  return(pval_list_sim)
  
}


pairedPowerAnalysis = function(params, sims=5, nmin, nmax, interval = 1) {
  
  result_matrix <- matrix(ncol = sims, nrow = 0)
  
  n = max(2, nmin)
  
  while (n <= nmax) {
    
    results = vector()
    
    for (randomseed in 1:sims) {
      
      simData = simulatePaired(params, n)
      
      pvalList = evalSimulatedData(simData, n)
      
      
      n_pos <- length(which(params$de & !is.na(pvalList[[1]])))
      
      n_tp <- length(which(pvalList[[1]] < 0.05 & params$de))
      
      results[randomseed] <- n_tp / n_pos
      
    }
    
    result_matrix <- rbind(result_matrix, results)
    
    rownames(result_matrix)[dim(result_matrix)[1]] <- n
    
    n <- n + interval
    
  }
  
  colnames(result_matrix) <- 1:sims
  
  return(result_matrix)
  
}





Data_wd <- "./Data/"


ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t", header=TRUE)[,-c(13)]

ko_noz_filt <- ko_noz[rowSums(ko_noz) > 12,]

ko_noz_filt <- filter(ko_noz_filt , rowSums(across(everything(), ~.x==0))<=9)


snf2_metadata <- read.table(paste(Data_wd, "metadata (1).txt", sep=""), sep="\t", fill = TRUE)

panelist = snf2_metadata$panelist

treatment = snf2_metadata$toothpaste



params = estimatePairedParams(ko_noz_filt, panelist, treatment)

results = pairedPowerAnalysis(params, 3, 5, 30, 5)


plot(rownames(results),rowMeans(results, na.rm=T), main="edgeR Simulation on Snf2 Data", xlab = "Number of Replicates", ylab = "Power")


