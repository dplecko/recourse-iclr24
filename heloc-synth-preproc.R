
# get the top6 variables + outcome
top6 <- c("ExternalRiskEstimate", "NetFractionRevolvingBurden", 
          "AverageMInFile", "PercentTradesWBalance",
          "PercentTradesNeverDelq", "MSinceMostRecentInqexcl7days")

out <- "RiskPerformance"

data <- read.csv("data/heloc_dataset_v1.csv")
data <- data[, c(top6, out)]
data[[out]] <- ifelse(data[[out]] == "Bad", 0, 1)

# remove the -9s -> how many samples are affected?
rm_idx <- which(rowMeans(do.call(cbind, lapply(data, function(x) x == -9))) > 0)
data <- data[-rm_idx, ]

# normalize things to percentages
trim_and_scale <- function(x) {
  
  x[x <= 0] <- 0 + runif(length(x[x <= 0]))
  x[x >= 100] <- 100 - runif(length(x[x >= 100]))
  x / 100
}

data$MSinceMostRecentInqexcl7days <- ifelse(
  data$MSinceMostRecentInqexcl7days < 0,
  0, data$MSinceMostRecentInqexcl7days
)

scale_vars <- c("ExternalRiskEstimate", "NetFractionRevolvingBurden",
                "PercentTradesWBalance", "PercentTradesNeverDelq")
for (var in scale_vars) {
  
  data[[var]] <- trim_and_scale(data[[var]])
}

write.csv(data, file = "data/heloc_top6.csv", row.names = FALSE)