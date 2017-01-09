library('coda')
library('data.table')
library('ggplot2')
library('grid')
library('gridExtra')
library('mvtnorm')
library('R2jags')
library('rjags')


PlotPriorPrNum <- function(K, Smax, hyper.pars.list, printout = TRUE) {
# Example:
#  PlotPriorPrNum(5, 5, hyper.pars.list)
  a = hyper.pars.list$pind_a
  b = hyper.pars.list$pind_b
  legend = paste0(c("E(theta1) = ", "E(theta2) = ", "SD(theta) = ",
                    "Pr1.alpha = ", "Pr1.beta = ",
                    "E(Pr.1) = ", "SD(Pr.1) = "),
                  c(hyper.pars.list$mu_logit[1], hyper.pars.list$theta2_mu,
                    round(sqrt(1/hyper.pars.list$tau_theta), 3),
                    a, b, a/(a + b), round(a * b/(a + b + 1)/(a + b) ^ 2, 3)))
  EtioPrior = ListEtiologyPriorSC1(K, Smax, hyper.pars.list,
                                   as.character(1:K), 20000)
  PrNum = as.data.frame(EtioPrior$Pr.num.pathogen)
  g = ggplot(data = PrNum) + xlab("Number of Pathogens") +
      ylab("Probability") + theme_bw() +
      geom_point(aes(x = 0:Smax, y = mean)) +
      geom_ribbon(aes(x = 0:Smax, ymin = 0, ymax = mean), alpha = 0.3) +
      geom_linerange(aes(x = 0:Smax, ymax = upper, ymin = lower),
                     alpha = 0.6) +
      geom_text(aes(x = (0:Smax) - 0.1, y = mean + 0.05,
                    label = as.character(mean)),
                angle = 90, vjust = 0) +
      annotate("text", x = Smax - 1.5, y = 0.85 - (1:length(legend)) * 0.05,
               hjust = 0, label = legend)
  if (printout) {
    print(g)
  }
  return(g)
}

PlotSimStudy <- function(all.fit, true.par) {
  # Example:
  #   true.par = TrueParVal(sim.study.2, par.to.save)
  #   PlotSimStudy(sim.study.2$stan.fit.all, true.par)
  dt = as.data.table(all.fit)
  dt = dt[!is.na(par.est)]
  dt = dt[!grepl("logit", par.name)]
  g = ggplot()
  print(
    g + geom_violin(data = dt, aes(x = method, y = par.est),
                    draw_quantiles = c(0.025, 0.5, 0.975)) +
      facet_grid(. ~ par.name) +
      geom_hline(data = true.par,
                 aes(yintercept = par.est, color = "red"))
  )
}


PlotCompareResult <- function(fit1, fit2, method.names = c("A", "B")) {
  # Example:
  #   PlotCompareResult(Mu.fit, eti_top5, c("SC1", "nlcm"))
  if (prod(colnames(fit1) == colnames(fit2)) < 1) {
    message("Variable names must be the same.")
    return(NULL)
  }
  dt = rbind(WideToLong(fit1, method.names[1]),
             WideToLong(fit2, method.names[2]))
  g = ggplot()
  print(
    g + geom_violin(data = dt, aes(x = Method, y = Estimate),
                    draw_quantiles = c(0.025, 0.5, 0.975)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_grid(. ~ Parameter) 
  )
}

MakeXLabel <- function(EtioMat, EtioComb, num.keep = NULL) {
  Label = apply(EtioMat, 1, function(x) {
    label = ""
    for (xi in x) {
      if (xi > 0) {
        label = paste(label, "X")
      } else {
        label = paste(label, "_")
      }
    }
    label
  })
  result  = data.frame(EtioComb = EtioComb, Label = Label)
  if (length(num.keep) > 0) {
    result = result[1:num.keep, ]
  }
  result
}

PlotByCombination <- function(coda.chains, sim.obj, hyper.pars.list,
                              etio.names = c("Patho_A", "Patho_B",
                                             "Patho_C", "Patho_D", "Patho_E"),
                              contrast = "prior", reorder = FALSE,
                              num.keep = NULL, baker.result = NULL,
                              has.true.value = FALSE) {
# Example:
# plog.obj = PlotByCombination(coda.fit[[1]], sim.obj, hyper.pars.list)
#  
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  EtioList = ListEtiology(coda.chains, sim.obj, etio.names, reorder, num.keep)
  EtioPrior = ListEtiologyPriorSC1(K, Smax, hyper.pars.list,
                                   etio.names, 20000, num.keep)
  X.unique = data.frame(uniquecombs(sim.obj$X))
  gplot.obj = list()
  for (i in 1:length(EtioList)) {
    Prior = merge(data.frame(EtioComb = EtioList[[i]]$EtioComb),
                  EtioPrior$Cell.prob, sort = FALSE)
    Label = MakeXLabel(EtioPrior$EtioMat, EtioPrior$Cell.prob$EtioComb,
                       num.keep)
    XLabel = merge(data.frame(EtioComb = EtioList[[i]]$EtioComb),
                   Label, sort = FALSE)
    legend  = paste0(c("Age = ", "Severity = "),
                     c(X.unique$AGE[i], X.unique$SEVERITY[i]))
    g = ggplot() + xlab("") +
        geom_point(data = EtioList[[i]], aes(x = 1:nrow(EtioList[[i]]),
                                             y = Probability)) +
        geom_linerange(data = EtioList[[i]], aes(x = 1:nrow(EtioList[[i]]),
                                                 ymax = Prob.upper,
                                                 ymin = Prob.lower)) +
        scale_x_continuous(breaks = 1:nrow(XLabel),
                           labels = as.character(XLabel$Label)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        annotate("text", x = nrow(EtioList[[i]]) - 4,
                 y = 0.7 - (1:length(legend)) * 0.1,
                 hjust = 0, label = legend)
        
    
    if (contrast == "prior") {
      g = g + geom_point(data = Prior, aes(x = (1:nrow(Prior)) + 0.2,
                                           y = Probability),
                         alpha = 0.5) +
              geom_linerange(data = Prior, aes(x = (1:nrow(Prior)) + 0.2,
                                               ymax = Prob.upper,
                                               ymin = Prob.lower),
                             alpha = 0.5)
    }
    if (contrast == "baker" && length(baker.result) > 0) {
      g = g + geom_point(data = baker.result[[i]],
                         aes(x = (1:nrow(baker.result[[i]])) + 0.2,
                             y = Probability),
                         alpha = 0.5) +
        geom_linerange(data = baker.result[[i]],
                       aes(x = (1:nrow(baker.result[[i]])) + 0.2,
                           ymax = Prob.upper, ymin = Prob.lower),
                       alpha = 0.5)      
    }
    if (has.true.value) {
      par0 = data.frame(Prob = sim.obj$cell.prob.unique[i, 1:num.keep])
      g = g + geom_point(data = par0, aes(x = 1:num.keep,
                                          y = Prob),
                         col = "red")
    }
    
    footnote <- paste(rev(etio.names), collapse = "\n")
    grid.newpage()
    gg = arrangeGrob(g, bottom = textGrob(
      footnote, x = 0.05, hjust = 0, vjust = -1.4,
      gp = gpar(fontface = "italic", fontsize = 7, lineheight = 0.95)))
    gplot.obj[[i]] = gg
    grid.draw(gg)
  }
  return(gplot.obj)
}

PlotByPathogen <- function(coda.chains, sim.obj,
                           etio.names = c("Patho_A", "Patho_B",
                                          "Patho_C", "Patho_D", "Patho_E"),
                           mu.fit = NULL) {
# Example:
# PlotByPathogen(coda.fit[[1]], sim.obj)
  Mu0 = sim.obj$pars.baseline$Mu
  Mu.true = data.frame(Value = as.vector(t(Mu0)),
                       Parameter = rep(etio.names, times = nrow(Mu0)),
                       Strata = rep(paste("strata", 1:nrow(Mu0)),
                                    each = ncol(Mu0)))
  if (length(mu.fit) > 0) {
    Mu.fit = mu.fit
  } else {
    Mu.fit = ExtractMu(coda.chains, sim.obj)
  }
  mu.samples = NULL
  for (i in 1:length(Mu.fit)) {
    mu.sample = Mu.fit[[i]]
    colnames(mu.sample) = etio.names
    mu.samples = rbind(mu.samples,
                       WideToLong(mu.sample, paste("strata", i)))
  }
  colnames(mu.samples)[2] = "Strata"
  g = ggplot()
  print(
    g + geom_violin(data = mu.samples, aes(x = Strata, y = Estimate),
                    draw_quantiles = c(0.025, 0.5, 0.975)) +
      facet_grid(. ~ Parameter) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") +
      geom_point(data = Mu.true,
                aes(x = Strata, y = Value),
                shape = 10,
                colour = "red")
  )
}

PlotPriorEffect <- function(sim.fit1, sim.fit2) {
  # TODO:
  Mu.fit1 = ExtractMuCellProb(sim.fit1)$Mu.fit
  Mu.fit2 = ExtractMuCellProb(sim.fit2)$Mu.fit
  
  K = ncol(Mu.fit1[[1]])
  
}