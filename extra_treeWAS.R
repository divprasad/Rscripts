

for (i in 1:length(sig.list)) {
  corr.dat <- sig.list[[i]]$corr.dat
  corr.sim <- sig.list[[i]]$corr.sim
  p.vals <- sig.list[[i]]$p.vals
  names(p.vals) <- names(corr.dat)
  sig.snps.names <- sig.list[[i]]$sig.snps.names
  sig.snps <- sig.list[[i]]$sig.snps
  sig.corrs <- sig.list[[i]]$sig.corrs
  sig.p.vals <- sig.list[[i]]$sig.p.vals
  min.p <- sig.list[[i]]$min.p
  sig.thresh <- sig.list[[i]]$sig.thresh
  if (plot.manhattan == TRUE) {
    manhattan.plot(p.vals = abs(corr.dat), col = "funky", 
                   transp = 0.25, sig.thresh = sig.thresh, thresh.col = "red", 
                   snps.assoc = NULL, snps.assoc.col = "red", jitter.amount = 1e-05, 
                   min.p = NULL, log10 = FALSE, ylab = paste(TEST[[i]], 
                                                             "score", sep = " "))
    title(paste("\n \n(", TEST[[i]], "score)"), cex.main = 0.9)
  }
  if (plot.null.dist == TRUE) {
    plot_sig_snps(corr.dat = abs(corr.dat), corr.sim = abs(corr.sim), 
                  corr.sim.subset = NULL, sig.corrs = abs(corr.dat[sig.snps]), 
                  sig.snps = sig.snps.names, sig.thresh = sig.thresh, 
                  test = TEST[[i]], sig.snps.col = "black", hist.col = rgb(0, 
                                                                           0, 1, 0.5), hist.subset.col = rgb(1, 0, 0, 
                                                                                                             0.5), thresh.col = "red", snps.assoc = NULL, 
                  snps.assoc.col = "blue", bg = "lightgrey", grid = TRUE, 
                  freq = FALSE, plot.null.dist = TRUE, plot.dist = FALSE)
  }
  if (plot.dist == TRUE) {
    plot_sig_snps(corr.dat = abs(corr.dat), corr.sim = abs(corr.sim), 
                  corr.sim.subset = NULL, sig.corrs = abs(corr.dat[sig.snps]), 
                  sig.snps = sig.snps.names, sig.thresh = sig.thresh, 
                  test = TEST[[i]], sig.snps.col = "black", hist.col = rgb(0, 
                                                                           0, 1, 0.5), hist.subset.col = rgb(1, 0, 0, 
                                                                                                             0.5), thresh.col = "red", snps.assoc = NULL, 
                  snps.assoc.col = "blue", bg = "lightgrey", grid = TRUE, 
                  freq = FALSE, plot.null.dist = FALSE, plot.dist = TRUE)
  }
  phen.curr <- phen
  if (length(sig.snps) == 0) 
    sig.snps <- sig.corrs <- NULL
  if (length(sig.snps) > 0) {
    toKeep <- sig.snps
    snps.toKeep <- snps[, toKeep]
    G1P1 <- G0P0 <- G1P0 <- G0P1 <- NA
    levs <- unique(as.vector(unlist(phen)))
    levs <- levs[!is.na(levs)]
    n.levs <- length(levs)
    if (n.levs == 2) {
      noms <- names(phen)
      phen <- as.numeric(as.factor(phen))
      phen <- rescale(phen, to = c(0, 1))
      names(phen) <- noms
      if (length(toKeep) > 1) {
        G1P1 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps.toKeep[which(phen == 
                                                                                            1), e] == 1)))
        G0P0 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps.toKeep[which(phen == 
                                                                                            0), e] == 0)))
        G1P0 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps.toKeep[which(phen == 
                                                                                            0), e] == 1)))
        G0P1 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps.toKeep[which(phen == 
                                                                                            1), e] == 0)))
      }
      else {
        G1P1 <- length(which(snps.toKeep[which(phen == 
                                                 1)] == 1))
        G0P0 <- length(which(snps.toKeep[which(phen == 
                                                 0)] == 0))
        G1P0 <- length(which(snps.toKeep[which(phen == 
                                                 0)] == 1))
        G0P1 <- length(which(snps.toKeep[which(phen == 
                                                 1)] == 0))
      }
      df <- data.frame(sig.snps, sig.p.vals, sig.corrs, 
                       G1P1, G0P0, G1P0, G0P1)
      names(df) <- c("SNP.locus", "p.value", "score", 
                     "G1P1", "G0P0", "G1P0", "G0P1")
    }
    else {
      df <- data.frame(sig.snps, sig.p.vals, sig.corrs)
      names(df) <- c("SNP.locus", "p.value", "score")
    }
  }
  else {
    df <- "No significant SNPs found."
  }
  min.p <- 1/length(corr.sim)
  names(min.p) <- c("p-values listed as 0 are less than:")
  results <- list()
  results[[1]] <- corr.dat
  results[[2]] <- corr.sim
  results[[3]] <- p.vals
  results[[4]] <- sig.thresh
  results[[5]] <- df
  results[[6]] <- min.p
  names(results) <- c("corr.dat", "corr.sim", "p.vals", 
                      "sig.thresh", "sig.snps", "min.p.value")
  RES[[i]] <- results
}


