library(MASS)
library(edgeR)
library(dplyr)

estimateOFParams = function(counts, treatment) {
  
  # sort data by treatment
  
  sort_order <- order(treatment)
  
  treatment <- treatment[sort_order]
  
  counts <- counts[, sort_order]
  
  # create the DGEList object using the sorted counts
  
  y = DGEList(counts=counts)
  
  y <- calcNormFactors(y)
  
  # create design matrix
  
  design = model.matrix(~treatment)
  
  rownames(design) <- colnames(y)
  
  # estimate dispersions
  
  dispsCR = dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, 
                                          trend = FALSE, AveLogCPM = NULL, min.row.sum = 6, prior.df = 0, span = 0.3, 
                                          grid.npts = 15, grid.range = c(-8, 8))
  
  sample_data <- data.frame(treatment)
  
  sample_data$libsize <- log(colSums(y$counts))
  
  libsize <- sample_data$libsize
  
  nofit <- c()
  
  fc <- matrix(nrow = dim(y$counts)[1], ncol = 2)
  
  for (i in 1:dim(y$counts)[1]) {
    f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
    tryCatch({ glm(y$counts[i,] ~ treatment + 0, offset = libsize, family = f) -> fit },
             warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
    fc[i,] <- fit$coefficients
  }
  
  y <- DGEList(counts = ko_noz_filt[-nofit,])
  
  params = list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
  
  
  # find de for given count data
  
  edgeR_cds <- DGEList(counts = counts[-nofit,], group = treatment)
  
  edgeR_cds <- calcNormFactors(edgeR_cds)
  
  edgeR_cds <- estimateCommonDisp(edgeR_cds)
  
  edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
  
  res <- exactTest(edgeR_cds, pair = c(unique(treatment)[1], unique(treatment)[2]))$table
  
  pval <- res$PValue
  
  de = pval < 0.05
  
  params[["de"]] <- de
  
  return(params)
  
}


simulateOF = function(params, n) {
  
  fcCopy = params$fc
  
  de = params$de
  
  libsize = params$sample_data$libsize
  
  # mean_expr <- (fcCopy[, 1] + fcCopy[, 2]) / 2
  
  # fcCopy[!de,] <- c(mean_expr[!de],mean_expr[!de])
  
  
  
  sim_libsizes <- runif(n * 2, min = min(libsize), max = max(libsize))
  
  
  m <- matrix(nrow = length(params$dispsCR), ncol = n * 2)

  
  for (i in 1:length(params$dispsCR)) {
    for (j in 1:(n * 2)) {
      m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] + ifelse(j <= n, fcCopy[i, 1], fcCopy[i, 2])), size = 1 / params$dispsCR[i])
    }
  }
  
  
  
  label <- paste0(c(rep("A", n), rep("B", n)), rep(1:n, 2))
  
  colnames(m) <- label
  
  rownames(m) <- paste0("g", 1:length(params$dispsCR))
  
  return(m)
  
}





evalOFData = function(m, n) {
  
  tryCatch({
    
    treatmentFactor <- factor(rep(c("N", "T"), n))
    
    edgeR_cds <- DGEList(m, group=treatmentFactor)
    
    edgeR_cds <- calcNormFactors(edgeR_cds)
    
    edgeR_cds <- estimateCommonDisp(edgeR_cds)
    
    edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
    
    res <- exactTest(edgeR_cds, pair = c(unique(treatmentFactor)[1], unique(treatmentFactor)[2]))$table
    
    pvalSim <- res$PValue
    
    padjSim <- p.adjust(pvalSim, method = "BH")
    
    resSim <- cbind(pvalSim, padjSim)
    
    pval_list_sim[["er"]] <- as.matrix(resSim)
    
  }, error = function(e) { printerror(e, "edgeR") })
  
  return(pval_list_sim)
  
}


ofPowerAnalysis = function(params, sims=5, nmin, nmax, interval = 1) {
  
  result_matrix <- matrix(ncol = sims, nrow = 0)
  
  n = max(2, nmin)
  
  while (n <= nmax) {
    
    results = vector()
    
    for (randomseed in 1:sims) {
      
      simData = simulateOF(params, n)
      
      pvalList = evalOFData(simData, n)
      
      
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

# unused in this one factor experiment

# panelist = snf2_metadata$panelist

treatment = snf2_metadata$toothpaste



paramsOF = estimateOFParams(ko_noz_filt, treatment)

resultsOF = ofPowerAnalysis(paramsOF, 3, 5, 10, 1)


plot(rownames(resultsOF),rowMeans(resultsOF, na.rm=T), main="edgeR Simulation on Snf2 Data", xlab = "Number of Replicates", ylab = "Power")


