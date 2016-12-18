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

MakeXLabel <- function(EtioMat, EtioComb) {
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
  data.frame(EtioComb = EtioComb, Label = Label)
}

PlotByCombination <- function(coda.chains, sim.obj, hyper.pars.list,
                              etio.names = c("Patho_A", "Patho_B",
                                             "Patho_C", "Patho_D", "Patho_E"),
                              contrast = "prior", reorder = FALSE) {
# Example:
# PlotByCombination(coda.fit[[1]], sim.obj, hyper.pars.list)
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  EtioList = ListEtiology(coda.chains, sim.obj, etio.names, reorder)
  EtioPrior = ListEtiologyPriorSC1(K, Smax, hyper.pars.list,
                                   etio.names, 20000)

  for (i in 1:length(EtioList)) {
    Prior = merge(data.frame(EtioComb = EtioList[[i]]$EtioComb),
                  EtioPrior$Cell.prob, sort = FALSE)
    Label = MakeXLabel(EtioPrior$EtioMat, EtioPrior$Cell.prob$EtioComb)
    XLabel = merge(data.frame(EtioComb = EtioList[[i]]$EtioComb),
                   Label, sort = FALSE)
    g = ggplot() + xlab("") +
        geom_point(data = EtioList[[i]], aes(x = 1:nrow(EtioList[[i]]),
                                             y = Probability)) +
        geom_linerange(data = EtioList[[i]], aes(x = 1:nrow(EtioList[[i]]),
                                                 ymax = Prob.upper,
                                                 ymin = Prob.lower)) +
        scale_x_continuous(breaks = 1:32,
                           labels = as.character(XLabel$Label)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
    
    if (contrast == "prior") {
      g = g + geom_point(data = Prior, aes(x = (1:nrow(Prior)) + 0.2,
                                           y = Probability),
                         alpha = 0.5) +
              geom_linerange(data = Prior, aes(x = (1:nrow(Prior)) + 0.2,
                                               ymax = Prob.upper,
                                               ymin = Prob.lower),
                             alpha = 0.5)
    }
    footnote <- paste(etio.names, collapse = "\n")
    grid.newpage()
    gg = arrangeGrob(g, bottom = textGrob(
      footnote, x = 0.05, hjust = 0, vjust = -1.4,
      gp = gpar(fontface = "italic", fontsize = 7, lineheight = 0.95)))
    grid.draw(gg)

  }
}

PlotByPathogen <- function(coda.chains, sim.obj,
                           etio.names = c("Patho_A", "Patho_B",
                                          "Patho_C", "Patho_D", "Patho_E")) {
# Example:
# PlotByPathogen(coda.fit[[1]], sim.obj)
  Mu0 = sim.obj$pars.baseline$Mu
  Mu.true = data.frame(Value = as.vector(t(Mu0)),
                       Parameter = rep(etio.names, times = nrow(Mu0)),
                       Strata = rep(paste("strata", 1:nrow(Mu0)),
                                    each = ncol(Mu0)))
  Mu.fit = ExtractMu(coda.chains, sim.obj)
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

