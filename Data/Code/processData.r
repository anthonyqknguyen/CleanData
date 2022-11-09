# Setting woking directory
# assumes starting location in terminal is /CleanData
setwd(getwd())


#install.packages("tidyr")
#install.packages("rlang")
library(tidyr)
sessionInfo()
library(rlang)
# Setting the data directory
Data_wd <- "./Data/"

# Reading in the counts table
counts <- read.table(paste(Data_wd, "counts.txt", sep=""), sep="\t", header =TRUE,
                     row.names = 1)


# Reading in the metadata
metadata <- read.delim(paste(Data_wd, "metadata.txt", sep =""), header = TRUE, sep = "\t", row.names = 1)
metadata_filt <- metadata[complete.cases(metadata[ , "Smoker"]),]

#counts_filt <- t(counts)[rownames(t(counts)) %in% rownames(metadata_filt)]

counts <- as.data.frame(t(counts))
#counts_filt <- counts[rownames(counts) %in% rownames(metadata_filt)]

counts_filt <- dplyr::select(counts, rownames(metadata_filt))
counts_filt <- counts_filt * 5
counts_filt[counts_filt == 0] <- 1


# Reading in the taxa data
taxa <- read.table(paste(Data_wd, "taxa.txt", sep=""), sep="\t")

# Reading in Khier data
khier <- read.table(paste(Data_wd, "khier_2020_fixed.txt", sep=""), sep="\t")

# Reading in ko_noz data
ko_noz <- read.table(paste(Data_wd, "ko_noz_rn_tmm.txt", sep=""), sep="\t")

# Reading in snf2 metadata
snf2_metadata <- read.table(paste(Data_wd, "metadata (1).txt", sep=""), sep="\t", fill = TRUE)

#View(taxa)
#View(metadata)
library("DESeq2")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = metadata_filt, 
                               design = ~Smoker)
res <- DESeq(dds, fitType = 'glmGamPoi')
install.packages("Matrix")
#View(dds$Smoker)

of_estimate_params <- function(rawdata, condition) {
  y <- DGEList(counts = rawdata)
  y <- calcNormFactors(y)
  
  design <- model.matrix(~condition)
  rownames(design) <- colnames(y)
  
  dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, trend = FALSE, AveLogCPM = NULL, min.row.sum = 5, prior.df = 0, span = 0.3, grid.npts = 15, grid.range = c(-8, 8)) -> dispsCR
  
  sample_data <- data.frame(condition)
  #sample_data$libsize = log(colSums(y$counts)) - log(nrow(y$counts))
  sample_data$libsize <- log(colSums(y$counts))
  libsize <- sample_data$libsize
  nofit <- 1000000
  fc <- matrix(nrow = dim(y$counts)[1], ncol = 2)
  for (i in 1:dim(y$counts)[1]) {
    f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
    tryCatch({ glm(y$counts[i,] ~ condition + 0, offset = libsize, family = f) -> fit },
             warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
    fc[i,] <- fit$coefficients
  }
  y <- DGEList(counts = rawdata[-nofit,])
  list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
}

of_DE_call <- function(rawdata, condition) {
  #DESeq2#
  dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(condition), design = ~condition)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  pval <- res$pval
  padj <- res$padj
  res <- cbind(pval, padj)
  ds2 <- as.matrix(res)
  rm(res, pval, padj)
  
  #packages = c("ds2", "ds","er","ss","eb")
  packages <- c("ds2", "er", "ss", "eb")
  
  de <- rep(TRUE, dim(rawdata)[1])
  for (i in packages) {
    temp <- length(which(get(i)[, "padj"] < 0.05))
    print(paste(i, ": number of DE called", temp))
    de <- de & get(i)[, "padj"] < 0.05
  }
  print(paste("intersection :", length(which(de))))
  de[is.na(de)] <- FALSE
  de
}

print_params <- function(params, de) {
  y <- params$y
  cat("number of genes:\t")
  cat(dim(y$counts)[1])
  cat("\n")
  temp <- cpm(y, log = TRUE, prior.count = 1)
  cat("cpm")
  cat("\t")
  cat(signif(median(temp), digits = 3))
  cat(" (")
  cat(signif(quantile(temp, 0.25), digits = 3))
  cat(" - ")
  cat(signif(quantile(temp, 0.75), digits = 3))
  cat(")\n")
  
  
  temp <- params$dispsCR
  cat("dispersion\t")
  cat(signif(median(temp), digits = 3))
  cat(" (")
  cat(signif(quantile(temp, 0.25), digits = 3))
  cat(" - ")
  cat(signif(quantile(temp, 0.75), digits = 3))
  cat(")\n")
  
  if (!is.null(params$fc)) {
    fc <- params$fc
    if (dim(fc)[2] == 2) {
      temp <- log2(exp(abs(fc[, 1] - fc[, 2])))
    } else {
      temp <- log2(exp(abs(fc[, dim(fc)[2]])))
    }
    temp <- temp[de]
    cat("fc\t")
    cat(signif(median(temp), digits = 3))
    cat(" (")
    cat(signif(quantile(temp, 0.25), digits = 3))
    cat(" - ")
    cat(signif(quantile(temp, 0.75), digits = 3))
    cat(")\n")
  }
  
  cat("libsize\t")
  if (is.null(params$sample_data)) {
    libsize <- params$libsize
  } else {
    libsize <- params$sample_data$libsize
  }
  libsize <- log10(exp(libsize))
  cat(signif(mean(libsize), digits = 3))
  cat(" +/- ")
  cat(signif(var(libsize), digits = 3))
  cat("\n")
}

of_simdata <- function(disps, libsizes, fc, n, de, randomseed = 1) {
  set.seed(randomseed)
  mean_expr <- (fc[, 1] + fc[, 2]) / 2
  #set fc truth data
  #fc[!de,] <- c(mean_expr[!de],mean_expr[!de])
  
  #generate library sizes for n samples based around the mean of the library sizes in the preliminary data
  sim_libsizes <- runif(n * 2, min = min(libsizes), max = max(libsizes))
  
  m <- matrix(nrow = length(disps), ncol = n * 2)
  
  for (i in 1:length(disps)) {
    for (j in 1:(n * 2)) {
      m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] + ifelse(j <= n, fc[i, 1], fc[i, 2])), size = 1 / disps[i])
    }
  }
  #m[m >= 10000] <- 10000
  label <- paste0(c(rep("A", n), rep("B", n)), rep(1:n, 2))
  colnames(m) <- label
  rownames(m) <- paste0("g", 1:length(disps))
  m
}

of_eval <- function(rawdata, condition, program) {
  pval_list <- list()
  
  if (program == "DESeq2") {
    #DESeq2#
    tryCatch({
      dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(condition), design = ~condition)
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds)
      res <- results(dds)
      pval <- res$pval
      padj <- res$padj
      res <- cbind(pval, padj)
      pval_list[["ds2"]] <- as.matrix(res)
      rm(dds, res, pval, padj)
      gc()
    }, error = function(e) { printerror(e, "DESeq2") })
  }
  #     else if (program == "DESeq") {
  # 	#DESeq#
  # 	tryCatch({
  # 		DESeq_cds <- newCountDataSet(rawdata, condition)
  # 		DESeq_cds <- estimateSizeFactors(DESeq_cds)
  # 		DESeq_cds <- estimateDispersions(DESeq_cds)
  # 		pval <- nbinomTest(DESeq_cds, "A", "B", pvals_only = TRUE)
  # 		padj <- p.adjust(pval, method = "BH")
  # 		res <- cbind(pval, padj)
  # 		pval_list[["ds"]] <- as.matrix(res)
  # 		rm(DESeq_cds, res, pval, padj)
  # 		gc()
  # 	}, error = function(e) { printerror(e, "DESeq") })
  # }
  else if (program == "edgeR") {
    #edgeR#
    tryCatch({
      edgeR_cds <- DGEList(rawdata, group = condition)
      edgeR_cds <- calcNormFactors(edgeR_cds)
      edgeR_cds <- estimateCommonDisp(edgeR_cds)
      edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
      res <- exactTest(edgeR_cds, pair = c("A", "B"))$table
      pval <- res$PValue
      padj <- p.adjust(pval, method = "BH")
      res <- cbind(pval, padj)
      pval_list[["er"]] <- as.matrix(res)
      rm(edgeR_cds, res, pval, padj)
      gc()
    }, error = function(e) { printerror(e, "edgeR") })
  } else if (program == "sSeq") {
    #sSeq#
    tryCatch({
      as.character(condition) -> sSeq_condition
      res <- nbTestSH(rawdata, sSeq_condition, condA = unique(sSeq_condition)[1], condB = unique(sSeq_condition)[2])
      pval <- res$pval
      padj <- p.adjust(pval, method = "BH")
      res <- cbind(pval, padj)
      pval_list[["ss"]] <- as.matrix(res)
      rm(res, sSeq_condition, pval, padj)
      gc()
    }, error = function(e) { printerror(e, "sSeq") })
  } else if (program == "EBSeq") {
    #EBSeq
    tryCatch({
      Sizes <- MedianNorm(rawdata)
      EBOut <- EBTest(Data = rawdata, Conditions = condition, sizeFactors = Sizes, maxround = 5)
      data.frame(pval = 1 - GetPP(EBOut)) -> temp0
      temp1 <- rawdata
      merge(temp1, temp0, all.x = TRUE, by.x = 0, by.y = 0) -> temp2
      pval <- temp2[, "pval"]
      names(pval) <- temp2[, "Row.names"]
      pval <- pval[rownames(rawdata)]
      padj <- pval
      res <- cbind(pval, padj)
      pval_list[["eb"]] <- as.matrix(res)
      rm(temp0, temp1, temp2, EBOut, Sizes, res, pval, padj)
      gc()
    }, error = function(e) { printerror(e, "EBSeq") })
  } else { stop("please select a program: DESeq2, edgeR, EBSeq or sSeq") }
  pval_list
}


RS_simulation <- function(budget = 3000, per_sample_price = 241, lane_size = 150e6, lane_price = 1331, mapping_proportion = 0.2,
                          sims = 5, params, designtype = "one factor", nmax = 20, nmin = 2, program) {
  attach(params)
  result_matrix <- matrix(ncol = sims, nrow = 0)
  n <- max(2, nmin)
  meanlibsize <- (lane_size * mapping_proportion) * (budget - (2 * n * per_sample_price)) / lane_price;
  while (meanlibsize > 100000 & n <= nmax) {
    results <- vector()
    for (randomseed in 1:sims) {
      #if meanlibsize specified (not zero) then change the mean libsize in the sample_data
      if (meanlibsize != "0") {
        temp <- log(meanlibsize)
        sample_data$libsize <- sample_data$libsize - (mean(max(sample_data$libsize), min(sample_data$libsize)) - temp)
      }
      if (designtype == "one factor") {
        m_sim <- of_simdata(disps = dispsCR, libsizes = sample_data$libsize, fc = fc, n = n, de, randomseed = randomseed)
        condition <- c(rep("A", n), rep("B", n))
        suppressMessages(pval_list <- of_eval(m_sim, condition, program))
      } else if (designtype == "paired") {
        m_sim <- paired_simdata(disps = dispsCR, libsizes = sample_data$libsize, fc_patient = fc[, 1:(dim(fc)[2] - 1)], fc = fc[, dim(fc)[2]], n = n, de, randomseed = randomseed)
        Patient <- factor(sort(c(1:n, 1:n)))
        Tissue <- factor(rep(c("N", "T"), n))
        suppressMessages(pval_list <- paired_eval(m_sim, Patient, Tissue))
      }
      n_pos <- length(which(de & !is.na(pval_list[[1]])))
      n_tp <- length(which(pval_list[[1]] < 0.05 & de))
      results[randomseed] <- n_tp / n_pos
    }
    result_matrix <- rbind(result_matrix, results)
    rownames(result_matrix)[dim(result_matrix)[1]] <- n
    n <- n + 1
    meanlibsize <- (lane_size * mapping_proportion) * (budget - (2 * n * per_sample_price)) / lane_price;
  }
  detach(params)
  colnames(result_matrix) <- 1:sims
  result_matrix
}

# only use of from github
# patient = metadata%panelist
# treatment = metadata$treatment
View(res)
