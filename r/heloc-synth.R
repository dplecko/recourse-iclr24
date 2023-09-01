
heloc_synth <- function(nsynth = 10^4, frank = TRUE, theta = 10,
                        last_n = 50, mc_samp = 1) {
  
  data <- read.csv(file.path(root, "data/heloc_top6.csv"))
  
  regr <- list(
    AverageMInFile = NULL, 
    PercentTradesWBalance = "beta", 
    PercentTradesNeverDelq = "beta", 
    NetFractionRevolvingBurden = "beta", 
    MSinceMostRecentInqexcl7days = "geom", 
    ExternalRiskEstimate = "beta"
  )
  
  top.ord <- names(regr)
  scm <- list()
  
  # have a pass retaining the objects
  for (var in top.ord[-1]) {
    
    scm[[var]][["object"]] <- reg_fun(var, top.ord, data, regr[[var]])
  }
  
  # have a pass generating data
  syn_dat <- data.frame(
    AverageMInFile = rnorm(nsynth, mean(data$AverageMInFile),
                           sd(data$AverageMInFile))
  )
  
  if (mc_samp > 1) {
    
    for (i in seq_len(mc_samp - 1)) {
      
      syn_dat <- rbind(syn_dat, tail(syn_dat, n = last_n))
    }
  }
  
  for (var in top.ord[-1]) {
    
    quants <- runif(nsynth)
    quants <- c(quants, rep(tail(quants, n = last_n), mc_samp - 1))
    scm[[var]][["U.pre"]] <- quants
    new_var <- gen_fun(scm[[var]][["object"]], syn_dat, quants)
    new_var <- data.frame(new_var)
    names(new_var) <- var
    syn_dat <- cbind(syn_dat, new_var)
  }
  
  # have a pass generating recourse data
  syn_rec <- syn_dat # can have some sort of interesting subsetting
  
  syn_rec$PercentTradesNeverDelq <- 0.85
  syn_rec$NetFractionRevolvingBurden <- 0.2

  for (var in c("MSinceMostRecentInqexcl7days", "ExternalRiskEstimate")) {
    
    syn_rec[[var]] <- NULL
    quants <- scm[[var]][["U.pre"]]
    
    if (frank) {
      
      quants <- BiCopCondSim(length(quants), cond.val = quants, cond.var = 1,
                             family = 5, par = theta)
    } else {
      
      # uniformly between current and 1!
      quants <- runif(length(quants), min = quants, max = 1)
    }
    
    scm[[var]][["U.post"]] <- quants
    new_var <- gen_fun(scm[[var]][["object"]], syn_rec, quants)
    new_var <- data.frame(new_var)
    names(new_var) <- var
    syn_rec <- cbind(syn_rec, new_var)
  }
  
  list(
    data = data,
    synth_data = syn_dat,
    rec_data = syn_rec,
    rec_mc = tail(syn_rec, n = last_n * mc_samp)
  )
}

