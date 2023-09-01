
library(VineCopula)

plt <- list()
for (tau in c(-0.5, 0.001, 0.5, 1)) {
  
  if (tau == 1) {
    
    x <- runif(200)
    df <- data.frame(X1 = x, X2 = x)
  } else {
    
    df <- data.frame(BiCopSim(500, family = 5, par = tau_to_theta(tau)))
  }
  
  plt[[as.character(tau)]] <- ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 0.2) + theme_minimal() +
    xlab(latex2exp::TeX("$Q_1$")) +
    ylab(latex2exp::TeX("$Q_2$")) +
    ggtitle(latex2exp::TeX(paste("$\\tau =", round(tau, 1), "$"))) +
    theme(legend.title = element_text(size = 14),
          axis.title = element_text(size = 12))
}

cowplot::plot_grid(plotlist = plt, labels = c("(a)", "(b)", "(c)", "(d)"),
                   ncol = 4L)

ggsave(file = "figures/frank_tau_dispersion.png", width = 10, height = 2)
