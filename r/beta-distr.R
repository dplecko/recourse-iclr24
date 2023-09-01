
expit <- function(x) exp(x) / (1 + exp(x))

beta_mom2par <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

geom_mom2par <- function(mu) 1 / (1 + mu)

tau_to_theta <- function(tau) BiCopTau2Par(5, tau)

theta_to_tau <- function(theta) BiCopPar2Tau(5, theta)
