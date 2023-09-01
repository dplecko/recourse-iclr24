library(mvtnorm)
library(ggplot2)

# set sample sizes
n0 <- n <- 500
nt <- 200

gen_Sigma <- function(p) {
  
  A <- matrix(rnorm(p * p), p, p)
  t(A) %*% A
}

p0 <- p <- 6

Sigma0 <- gen_Sigma(p0)
Sigma <- gen_Sigma(p)

beta0 <- runif(p0, -1, 1)
beta <- runif(p, -1, 1)

df <- pef <- NULL
nrep <- 500
for (ro in seq(0.1, 0.9, 0.1)) {
  
  for (i in seq_len(nrep)) {
    
    X0 <- rmvnorm(n0 + nt, sigma = Sigma0)
    X <- rmvnorm(n + nt, sigma = Sigma)
    
    # generate the noise
    epsAB <- rmvnorm(n+nt, sigma = matrix(c(1, ro, ro, 1), ncol = 2))
    
    # stage I: Y_0 = X_0 b_0 + eps_0
    Y0 <- X0 %*% beta0 + epsAB[, 1]
    eps_hat <- Y0 - lm(Y0 ~ X0)$fitted.values
    
    # stage II: Y = Xb + eps
    Y <- X %*% beta + epsAB[, 2]
    mod_e <- lm(Y ~ ., data.frame(Y = Y[seq_len(n)], X = X[seq_len(n), ], 
                                  eps_hat = eps_hat[seq_len(n)]))
    
    mod <- lm(Y ~ ., data.frame(Y = Y[seq_len(n)], X = X[seq_len(n), ]))
    
    df <- rbind(
      df, 
      data.frame(coef = names(coef(mod_e)[-1]), 
                 value = coef(mod_e)[-1], 
                 rep = i, rho = ro)
    )
    
    # MSPE
    Yhat_e <- predict(mod_e, data.frame(X = X[-seq_len(n), ], 
                                        eps_hat = eps_hat[-seq_len(n)]))
    
    Yhat <- predict(mod, data.frame(X = X[-seq_len(n), ]))
    
    pe <- mean((Y[-seq_len(n)] - Yhat)^2)
    pe_e <- mean((Y[-seq_len(n)] - Yhat_e)^2)
    pef <- rbind(pef, data.frame(PE = pe, PE_extended = pe_e, rep = i, rho = ro))
  }
}

# Fig. 1
ggplot(df[df$coef != "eps_hat",], aes(x = value)) +
  geom_density() + theme_minimal() +
  facet_wrap(~ coef, ncol = 3, nrow = 2, scales = "free_y") +
  geom_vline(
    data = data.frame(coef = sort(unique(df[df$coef != "eps_hat",]$coef)), 
                      value = c(beta)),
    aes(xintercept = value)
  ) +
  xlab(latex2exp::TeX("$\\widehat{\\beta}'_i$ Estimate")) +
  ylab("Density")

ggsave("figures/beta-hat.png", height = 3, width = 6 * 1.41)

# Fig. 2
ggplot(df[df$coef == "eps_hat",], 
       aes(x = value, fill = factor(rho))) +
  geom_density() + theme_minimal() +
  scale_fill_discrete(name = latex2exp::TeX("$\\rho$ value")) +
  # theme(legend.position = "bottom") +
  xlab(latex2exp::TeX("$\\widehat{\\beta}(\\widehat{\\epsilon})$ Estimate")) +
  ylab("Density")

ggsave("figures/beta-eps-hat.png", height = 3, width = 3 * 1.41)

# Fig. 3
pef$eff <- pef$PE_extended / pef$PE
ggplot(pef, aes(x = rho, y = eff, group = rho)) +
  geom_boxplot() +
  theme_minimal() +
  geom_point(aes(x = rho, y = 1 - rho^2), color = "red") +
  geom_line(data = data.frame(rho = seq(0.1, 0.9, 0.1), grp = 1),
            aes(x = rho, y = 1 - rho^2, group = grp), color = "red") +
  xlab(latex2exp::TeX("$\\rho$")) +
  ylab("MSPE(2-Stage) / MSPE(OLS)")

ggsave("figures/mspe-ratio.png", height = 3, width = 3 * 1.41)
