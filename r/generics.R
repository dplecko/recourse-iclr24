
visualizeGraph <- function(x, ...) {
  UseMethod("visualizeGraph", x)
}

visualizeGraph.fairadapt <- function(x, ...) plot(x$graph, ...)

predict.fairecourse <- function(object, data, rec_act, mc_samp = 1, 
                                use_obs = TRUE, ...) {

  rec_data <- data
  
  # if single recourse action, expand the list
  if (length(rec_act) == 1) {
    
    rec_act <- lapply(seq_len(nrow(data)), function(i) rec_act[[1]])
  } 
  
  # update the dataset based on recourse actions
  for (i in seq_along(rec_act)) {
    
    for (var in names(rec_act[[i]])) {
      
      rec_data[[var]] <- rec_act[[i]][[var]] 
    }
  }
  
  if (mc_samp > 1) {
    
    rec_data <- do.call(rbind, lapply(seq_len(mc_samp), function(x) rec_data))
    data <- do.call(rbind, lapply(seq_len(mc_samp), function(x) data))
  }
  
  engine <- object$q.engine
  
  for (var in setdiff(names(engine), all.vars(object$formula)[1L])) {
    
    # get the index of variables that need to be adjusted
    up.id <- is_updated(rec_act, var, object$adj.mat)
    
    assert_that(class(data[, var]) == engine[[var]][["type"]],
                msg = "Mismatch in column type with training data")

    # i) encode the variable if needed
    if (engine[[var]][["discrete"]]) {

      assert_that(
        all(data[, var] %in% engine[[var]][["unique.values"]]) |
          is.integer(engine[[var]][["discrete"]]),
        msg = paste0("New, unseen values of variable ", var, ". Disallowed.")
      )

      if (is.logical(engine[[var]][["discrete"]])) {
        
        data[, var] <- factor(data[, var],
                              levels = engine[[var]][["unique.values"]])
      }

      data[, var] <- as.integer(data[, var]) + 
                     runif(length(data[, var]), -0.5, 0.5)
    }

    # ii) computeQuants() based on data
    if (sum(up.id) > 0L) {
      
      if (!is.null(engine[[var]][["object.dual"]])) { # version C
        
        # create data.frame with _pre and _post columns
        post.dat <- rec_data[up.id,]
        pre.dat <- data[up.id,]
        names(post.dat) <- paste0(names(post.dat), "_post")
        names(pre.dat) <- paste0(names(pre.dat), "_pre")
        
        # generate new predictions
        if (use_obs) var.eng <- engine[[var]][["object.dual"]] else
          var.eng <- engine[[var]][["object.shadowdual"]]
        
        w.distr <- predict(var.eng, 
                           data = cbind(pre.dat, post.dat),
                           type = "quantiles", what = identity)$predictions
        
        # sample quantiles uniformly
        U.new <- runif(nrow(w.distr))
        w.hat <- vapply(seq_along(U.new), 
                        function(x) quantile(w.distr[x, ], U.new[x]),
                        numeric(1L))
        
        # assign the new value to recourse data.frame
        rec_data[up.id, var] <- w.hat
      } else if (!is.null(engine[[var]][["tau"]])) { # version B
        
        rec_data[up.id, var] <-
          computeQuants(
            engine[[var]][["object"]],
            data[, c(var, engine[[var]][["parents"]]), drop = FALSE],
            rec_data[up.id, engine[[var]][["parents"]], drop = FALSE],
            up.id, test = TRUE, tau = engine[[var]][["tau"]], ...
          )
      } else { # version A
        
        rec_data[up.id, var] <-
          computeQuants(
            engine[[var]][["object"]],
            data[, c(var, engine[[var]][["parents"]]), drop = FALSE],
            rec_data[up.id, engine[[var]][["parents"]], drop = FALSE],
            up.id, test = TRUE, ...
          )
      }
    }
    
    # iii) decode discrete
    if (engine[[var]][["discrete"]]) {

      if (is.integer(engine[[var]][["discrete"]])) {
        
        data[, var] <- as.integer(round(data[, var]))
        rec_data[, var] <- as.integer(round(rec_data[, var]))
      } else {
        
        data[, var] <-
          decodeDiscrete(data[, var], engine[[var]][["unique.values"]],
                         engine[[var]][["type"]], length(data[, var]))
        rec_data[, var] <-
          decodeDiscrete(rec_data[, var], engine[[var]][["unique.values"]],
                         engine[[var]][["type"]], length(rec_data[, var]))
      }
    }
  }

  if (mc_samp > 1) {
    
    rec_data <- cbind(rec_data, indiv = rep(seq_len(nrow(rec_data) / mc_samp), 
                                            mc_samp))
  }
  
  rec_data
}

