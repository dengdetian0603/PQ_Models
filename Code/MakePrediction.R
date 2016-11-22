library(Rcpp)


PredictCellProb <- function(X.new, MBS.new, MSS.new, patho.names,
							              K, Smax, X.train, coda.chains) {
	X.unique = uniquecombs(X.train)
	label = apply(X.unique, 1,
	              function(x) paste(as.character(x), collapse = ""))
	X.strata = which(label == paste(as.character(X.new), collapse = ""))

	design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Lall = rbind(rep(0, K), design.mat$Lmat)

  D = ncol(X.train)
  if (D == 1) {
    cell_prob = coda.chains[, paste0("cell_prob[", 1:(design.mat$J1 + 1), "]")]
  } else {
    cell_prob = coda.chains[, paste0("cell_prob[",
                                       X.strata, ",",
                                       1:(design.mat$J1 + 1), "]")]
  }
  A = apply(cell_prob, 1, sum)
  log_cell_prob = log(cell_prob/A)
  ss_tpr = coda.chains[, paste0("ss_tpr[", 1:K, "]")]
  bs_tpr = coda.chains[, paste0("bs_tpr[", 1:K, "]")]
  bs_fpr = coda.chains[, paste0("bs_fpr[", 1:K, "]")]

	pred.prob.dist = exp(predictLogProbMat(K, Lall, MSS.new, MBS.new,
                                         log_cell_prob,
                                         ss_tpr, bs_tpr, bs_fpr))
	pred.prob.dist = pred.prob.dist/rowSums(pred.prob.dist)
	
  pred.Mu.dist = pred.prob.dist[, -1] %*% t(design.mat$MuMat)
  pred.prob = apply(pred.prob.dist, 2, mean)
  pred.Mu = apply(pred.Mu.dist, 2, mean)

  list(pred.prob.dist = pred.prob.dist,
       pred.Mu.dist = pred.Mu.dist,
       pred.prob = pred.prob,
       pred.Mu = pred.Mu)
}