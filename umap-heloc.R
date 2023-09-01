
library(pracma)
library(umap)

# load heloc
data <- read.csv("data/heloc_dataset_v1.csv")
top10 <- c("ExternalRiskEstimate", "NetFractionRevolvingBurden", "AverageMInFile", 
           "MSinceOldestTradeOpen", "PercentTradesWBalance", "PercentInstallTrades", 
           "NumSatisfactoryTrades", "NumTotalTrades", "PercentTradesNeverDelq", 
           "MSinceMostRecentInqexcl7days")
out <- "RiskPerformance"
data <- data[, c(out, top10)]
data$RiskPerformance <- ifelse(data$RiskPerformance == "Bad", 0, 1)

# scale the columns to standard normals
for (i in 2:ncol(data)) {
  
  new_val <- data[, i]
  new_val <- (new_val - mean(new_val)) / sd(new_val)
  data[, i] <- new_val
}

# construct logistic model
logreg <- glm(RiskPerformance ~ ., data = data, family = "binomial")
probs <- predict(logreg, data, type = "response")

# plot the probabilities
ggplot(data.frame(x = probs), aes(x = x)) + geom_density() + theme_minimal()

# decide the decision boundary; handle the intercept; obtain normal vector
b.id <- which.min(abs(probs - 0.5))
x.b <- as.vector(as.matrix(data[b.id, -1]))


# pracma the normal vector
ort.cmp <- null(t(coef(logreg)[-1]))

# introduce perturbations of the point at boundary
n.bnd <- 3000
dirs <- do.call(cbind, lapply(seq_len(n.bnd), function(i) {
  prtb <- runif(9, -1, 1)
  runif(1, 0, 3) * prtb / sqrt(sum(prtb^2))
}))

bnd <- t(ort.cmp %*% dirs + x.b)

# umap
ump <- umap(as.matrix(data[, -1]))

# plot + color boundary (red, blue, green)
clr <- c(ifelse(probs < 0.5, -1, 1), rep(0, n.bnd))

ggplot(
  data.frame(rbind(ump$layout, predict(ump, bnd)), clr = clr),
  aes(x = X1, y = X2, color = factor(clr))
) +
  geom_point(size = ifelse(clr == 0, 0.5, 0.1)) + theme_minimal() +
  scale_color_manual(
    name = "Class", 
    values = c("#E41A1C", "#377EB8", "#4DAF4A"),
    labels = c("Negative", "Boundary", "Positive")) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 7))
