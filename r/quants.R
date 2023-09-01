
rangerQuants <- function(data, A.root, ind, min.node.size = 20, ...) {

  if (A.root) {
    return(
      structure(
        list(
          class0 = rangerQuants(data[ind, ], FALSE, NULL,
                                min.node.size = min.node.size, ...),
          class1 = rangerQuants(data[!ind, ], FALSE, NULL,
                                min.node.size = min.node.size, ...)),
        class = "rangersplit"
      )
    )
  }

  ranger::ranger(formula(data), data = data, quantreg = TRUE,
                 keep.inbag = TRUE, min.node.size = min.node.size, ...)
}


linearQuants <- function(data, A.root, ind,
                         tau = c(0.001, seq(0.005, 0.995, by = 0.01), 0.999),
                         ...) {

  if (A.root) {
    return(
      structure(
        list(
          class0 = linearQuants(data[ind, ], FALSE, NULL, tau = tau, ...),
          class1 = linearQuants(data[!ind, ], FALSE, NULL, tau = tau, ...)
        ),
        class = "quantregsplit"
      )
    )
  }

  offending.cols <- 1L + which(
    vapply(seq_col(data)[-1L], function(x) length(unique(data[, x])), 1L) == 1
  )

  keep.cols <- which(!(seq_col(data) %in% offending.cols))

  if (length(offending.cols) == (ncol(data) - 1L)) {
    form <- as.formula(paste(names(data)[1L], "~ 1"))
  } else {
    form <- formula(data[, keep.cols])
  }

  quantreg::rq(form, data = data, tau = tau, ...)
}

computeQuants <- function(x, data, newdata, ind, ...) {
  UseMethod("computeQuants", x)
}

computeQuants.ranger <- function(x, data, newdata, ind, test = FALSE, 
                                 emp.only = FALSE, tau = 1, ...) {

  # GetQuants
  if (isTRUE(test)) {

    empirical <- predict(x, data = data[, -1L, drop = FALSE],
                         type = "quantiles", what = identity)
    
    empirical <- empirical$predictions
    if (emp.only) return(empirical)

  } else {

    empirical <- x$random.node.values.oob
  }

  quantiles <- predict(x, data = newdata, type = "quantiles",
                       what = identity)
  quantiles <- quantiles$predictions
  
  inferQuant(data, empirical, quantiles, ind, tau)
}


#' @export
computeQuants.rqs <- function(x, data, newdata, ind, emp.only = FALSE, ...) {
  
  empirical <- predict(x, newdata = data[, -1L, drop = FALSE])
  if (emp.only) return(empirical)
  quantiles <- predict(x, newdata = newdata)

  inferQuant(data, empirical, quantiles, ind)
}

inferQuant <- function(data, empirical, quantiles, ind, tau) {

  eval <- data[, 1L]
  U.hat <- vapply(seq_row(data), function(x) ecdf(empirical[x, ]) (eval[x]),
                  numeric(1L))
  
  #' * new quantile is sampled from a copula. *
  if (is.null(tau) || tau == 1) {
    
    newU <- U.hat[ind]
  } else if (tau == -1) {
    
    newU <- 1 - U.hat[ind]
  } else {
    
    theta <- BiCopTau2Par(family = 5, tau = tau)
    newU <- BiCopCondSim(length(U.hat[ind]), cond.val = uwrap(U.hat[ind]), 
                         cond.var = 1, family = 5, par = theta)
  }
  
  vapply(seq_row(quantiles), function(x) quantile(quantiles[x, ], newU[x]),
         numeric(1L))
}
