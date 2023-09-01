library(ranger)

data <- read.csv("data/heloc_dataset_v1.csv")

head(data)

ncol(data)

data$RiskPerformance <- ifelse(data$RiskPerformance == "Bad", 0, 1)

out <- "RiskPerformance"
rf <- ranger(RiskPerformance ~ ., data = data, classification = TRUE,
             importance = "impurity")

top10 <- names(sort(importance(rf), decreasing = TRUE)[1:10])


rf_top <- ranger(RiskPerformance ~ ., data = data[, c(out, top10)], 
                 classification = TRUE, importance = "impurity")

rf_top

# run FCI on the data
reticulate::use_python("/usr/local/bin/python3")
library(reticulate)
py_config()$python

data_top <- r_to_py(data[, c(top10)])

# run python code

plot(fairadapt::graphModel(adj.mat))

# get the object

# --------------------- #

top10 <- c("ExternalRiskEstimate", "NetFractionRevolvingBurden", "AverageMInFile", 
           "MSinceOldestTradeOpen", "PercentTradesWBalance", "PercentInstallTrades", 
           "NumSatisfactoryTrades", "NumTotalTrades", "PercentTradesNeverDelq", 
           "MSinceMostRecentInqexcl7days")
out <- "RiskPerformance"
data[, c(top10)]

plt <- list()
for (var in top10) {
  
  plt[[var]] <- ggplot(data, aes_string(x = var)) +
    geom_density() + theme_minimal()
}
cowplot::plot_grid(plotlist = plt, ncol = 2)

coord <- c(1:20)
plot(fairadapt::graphModel(adj.mat), 
     layout = matrix(1:20, ncol = 2))

# --------------------- #
top6 <- c("ExternalRiskEstimate", "NetFractionRevolvingBurden", 
          "AverageMInFile", "PercentTradesWBalance",
          "PercentTradesNeverDelq", "MSinceMostRecentInqexcl7days")

adj.mat <- array(0, dim = c(7, 7))
colnames(adj.mat) <- rownames(adj.mat) <- c(top6, out)

grp <- list(
  g1 = c("AverageMInFile"),
  g3 = c("PercentTradesWBalance", "PercentTradesNeverDelq"),
  g4 = c("MSinceMostRecentInqexcl7days", "NetFractionRevolvingBurden"),
  g5 = c("ExternalRiskEstimate"),
  g6 = out
)

coord <- matrix(c(0, 0), ncol = 2L)

for (i in seq.int(2, length(grp))) {
  
  adj.mat[unlist(grp[seq.int(1, i-1)]), grp[[i]]] <- 1
  coord <- rbind(coord, cbind(c(seq_along(grp[[i]]) - mean(seq_along(grp[[i]]))), -(i-1)))
}

# sort topologically
top.ord <- topologicalOrdering(adj.mat)
adj.mat <- adj.mat[top.ord, top.ord]

plot(fairadapt::graphModel(adj.mat), layout = coord, vertex.size=20, 
     vertex.label.cex=0.8)

plt <- list()
for (var in top6) {
  
  plt[[var]] <- ggplot(data, aes_string(x = var)) +
    geom_density() + theme_minimal()
}
cowplot::plot_grid(plotlist = plt, ncol = 2)