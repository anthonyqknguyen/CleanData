library(MASS)
library(edgeR)
library(dplyr)

estimateOFParams = function(counts, factor) {
  
  # sort data by treatment
  
  sort_order <- order(factor)
  
  factor <- factor[sort_order]
  
  counts <- counts[, sort_order]
  
  # create the DGEList object using the sorted counts
  
  y = DGEList(counts=counts)
  
  y <- calcNormFactors(y)
  
  # create design matrix
  
  design = model.matrix(~factor)
  
  rownames(design) <- colnames(y)
  
  # estimate dispersions
  
  dispsCR = dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, 
                                          trend = FALSE, AveLogCPM = NULL, min.row.sum = 6, prior.df = 0, span = 0.3, 
                                          grid.npts = 15, grid.range = c(-8, 8))
  
  sample_data <- data.frame(factor)
  
  sample_data$libsize <- log(colSums(y$counts))
  
  libsize <- sample_data$libsize
  
  nofit <- c()
  
  fc <- matrix(nrow = dim(y$counts)[1], ncol = 2)
  
  for (i in 1:dim(y$counts)[1]) {
    f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
    tryCatch({ glm(y$counts[i,] ~ factor + 0, offset = libsize, family = f) -> fit },
             warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
    fc[i,] <- fit$coefficients
  }
  
  y <- DGEList(counts = counts[-nofit,])
  
  params = list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
  
  
  # find de for given count data
  
  edgeR_cds <- DGEList(counts = counts[-nofit,], group = factor)
  
  edgeR_cds <- calcNormFactors(edgeR_cds)
  
  edgeR_cds <- estimateCommonDisp(edgeR_cds)
  
  edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
  
  res <- exactTest(edgeR_cds, pair = c(unique(factor)[1], unique(factor)[2]))$table
  
  pval <- res$PValue
  
  de = pval < 0.05
  
  params[["de"]] <- de
  
  return(params)
  
}


estimatePairedParams = function(counts, factor1, factor2) {
  
  # sort data by treatment
  
  sort_order <- order(factor1)
  
  factor1 <- factor1[sort_order]
  
  counts <- counts[, sort_order]
  
  # create the DGEList object using the sorted counts
  
  y = DGEList(counts=counts)
  
  # create design matrix
  
  design = model.matrix(~factor2 + factor1)
  
  rownames(design) <- colnames(y)
  
  # estimate dispersions
  
  dispsCR = dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, 
                                          trend = FALSE, AveLogCPM = NULL, min.row.sum = 6, prior.df = 0, span = 0.3, 
                                          grid.npts = 15, grid.range = c(-8, 8))
  
  sample_data <- data.frame(factor2, factor1)
  
  sample_data$libsize <- log(colSums(y$counts))
  
  libsize <- sample_data$libsize
  
  nofit <- c()
  
  fc <- matrix(nrow = dim(y$counts)[1], ncol = length(levels(factor(factor2))) + length(levels(factor(factor1))) - 1)
  
  for (i in 1:dim(y$counts)[1]) {
    f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
    tryCatch({ glm(y$counts[i,] ~ factor2 + factor1 + 0, offset = libsize, family = f) -> fit },
             warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
    fc[i,] <- fit$coefficients
  }
  
  y <- DGEList(counts = counts[-nofit,])
  
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


simulateOF = function(params, n) {
  
  fc = params$fc
  
  de = params$de
  
  disps = params$dispsCR
  
  libsize = params$sample_data$libsize
  
  # mean_expr <- (fc[, 1] + fc[, 2]) / 2
  
  # fc[!de,] <- c(mean_expr[!de],mean_expr[!de])
  
  
  
  sim_libsizes <- runif(n * 2, min = min(libsize), max = max(libsize))
  
  
  m <- matrix(nrow = length(disps), ncol = n * 2)
  
  
  for (i in 1:length(disps)) {
    for (j in 1:(n * 2)) {
      m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] + ifelse(j <= n, fc[i, 1], fc[i, 2])), size = 1 / disps[i])
    }
  }
  
  m[m >= 50000] <- 50000
  
  label <- paste0(c(rep("A", n), rep("B", n)), rep(1:n, 2))
  
  colnames(m) <- label
  
  rownames(m) <- paste0("g", 1:length(disps))
  
  return(m)
  
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
  
  m[m >= 50000] <- 50000
  
  Panelist <- sort(c(1:n, 1:n))
  
  Treatment <- rep(c("N", "T"), n)
  
  colnames(m) <- paste0("s", Panelist, Treatment)
  
  rownames(m) <- paste0("g", 1:length(params$dispsCR))
  
  return(m)
  
}



evalOFData = function(m, n) {
  
  pval_list_sim = list()
  
  treatmentFactor <- factor(c(rep("N", n), rep("T", n)))
  
  edgeR_cds <- DGEList(m, group=treatmentFactor)
  
  edgeR_cds <- calcNormFactors(edgeR_cds)
  
  edgeR_cds <- estimateCommonDisp(edgeR_cds)
  
  edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
  
  res <- exactTest(edgeR_cds, pair = c(unique(treatmentFactor)[1], unique(treatmentFactor)[2]))$table
  
  pvalSim <- res$PValue
  
  padjSim <- p.adjust(pvalSim, method = "BH")
  
  resSim <- cbind(pvalSim, padjSim)
  
  pval_list_sim[["er"]] <- as.matrix(resSim)
  
  
  return(pval_list_sim)
  
}


evalPairedData = function(m, n) {
  
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


ofPowerAnalysis = function(counts, factor, sims = 5, nmin = 2, nmax = 10, interval = 1) {
  
  params = estimateOFParams(counts, factor)
  
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


pairedPowerAnalysis = function(counts, factor1, factor2, sims = 5, nmin = 2, nmax = 10, interval = 1) {
  
  params = estimatePairedParams(counts, factor1, factor2)
  
  result_matrix <- matrix(ncol = sims, nrow = 0)
  
  n = max(2, nmin)
  
  while (n <= nmax) {
    
    results = vector()
    
    for (randomseed in 1:sims) {
      
      simData = simulatePaired(params, n)
      
      pvalList = evalPairedData(simData, n)
      
      
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
