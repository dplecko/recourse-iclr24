
cnd_frank <- function(theta, u = 1/2, v = 1/2) {
  
  t1 <- exp(-theta * (u)) - 1
  t2 <- exp(-theta * (v)) - 1
  t3 <- exp(-theta) - 1
  
  -theta * exp(-theta * (u+v)) * t3 * (t3 + t1 * t2)^(-2)
}

infer_theta <- function(u, v) {
  
  th_seq <- seq(-35, 35, length.out = 500)
  log_l <- vapply(
    th_seq, function(theta) mean(log(cnd_frank(theta, u, v))), numeric(1L)
  )
  
  th_seq[which.max(log_l)]
} 

# ---------------- #
#  Copula Testing  #
# ---------------- #

frank_testing <- function(u, v, theta, browse = FALSE) {
  
  # if (browse) browser()
  # plot(u, v, pch = 19, cex = 0.2, xlim = c(0, 1), ylim = c(0, 1))
  # abline(0, 1, col = "red")
  nboot <- 500
  log_lt <- mean(log(cnd_frank(theta, u, v)))
  log_lcmp <- rep(0, nboot)
  for (i in seq_len(nboot)) {
    
    v_boot <- BiCopCondSim(length(u), cond.val = uwrap(u), cond.var = 1, 
                           family = 5, par = theta)
    th_boot <- infer_theta(u, v_boot)
    
    log_lcmp[i] <- mean(log(cnd_frank(th_boot, u, v_boot)))
  }
  
  mean(log_lcmp > log_lt)
}
