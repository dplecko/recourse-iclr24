
fairecourse <- function(formula, adj.mat, train.data, pre.rec = NULL,
                        post.rec = NULL, rec_act = NULL, cfd.mat = NULL, 
                        quant.method = rangerQuants,
                        visualize.graph = FALSE, eval.qfit = NULL, 
                        soft.test = FALSE, ...) {

  if (missing(adj.mat)) {
    adj.mat <- NULL
  }
  
  # if single recourse action, expand the list
  if (length(rec_act) == 1) {
    
    rec_act <- lapply(seq_len(nrow(pre.rec)), function(i) rec_act[[1]])
  } 
  
  # verify correctness of input
  # correctInput(formula, train.data, test.data, adj.mat, cfd.mat, top.ord,
  #              prot.attr, res.vars, quant.method)

  # reorder the adjacency matrix if necessary
  adj.mat <- adj.mat[colnames(adj.mat), ]

  # reorder columns and Factor-ize
  #' * check if formula argument can be dropped entirely *
  # train.data <- model.frame(formula, train.data)
  org.data <- train.data
  
  copula <- TRUE
  if (is.null(org.data)) {
    
    copula <- FALSE
    org.data <- pre.rec
  }

  # keep important parts of adjacency matrix
  adj.mat <- adj.mat[colnames(org.data), colnames(org.data)]

  if (is.null(cfd.mat) && !is.null(adj.mat)) {

    cfd.mat <- adj.mat
    cfd.mat[, ] <- 0
  }

  # construct the initial version of adapted data
  top.ord <- topologicalOrdering(adj.mat)
  
  ig <- graphModel(adj.mat, cfd.mat)

  q.engine <- list()

  # main procedure part
  var.ind <- seq.int(1L, length(top.ord))

  for (var in top.ord[var.ind]) {

    changed.parents <- getParents(var, adj.mat, top.ord)
    
    if (length(changed.parents) == 0) next

    q.engine[[var]] <- list()

    type <- class(org.data[, var])
    q.engine[[var]][["type"]] <- type

    discrete <- FALSE
    q.engine[[var]][["discrete"]] <- discrete

    curr.parents <- adjustmentSet(var, adj.mat, cfd.mat, top.ord)
    q.engine[[var]][["parents"]] <- curr.parents
    
    # check if Discrete
    if (length(unique(org.data[, var])) < 10 |
        is.factor(org.data[, var]) | is.integer(org.data[, var]) |
        is.character(org.data[, var])) {
      
      discrete <- TRUE
      q.engine[[var]][["discrete"]] <- discrete

      if (is.character(org.data[, var])) {
        org.data[, var] <- factor(org.data[, var])
      }

      if (is.factor(org.data[, var])) {
        
        org.data[, var] <- factor(
          org.data[, var], 
          levels = catOrder(org.data[, 1L], org.data[, var])
        )
        if (!is.null(pre.rec) && !is.null(post.rec)) {
          
          pre.rec[, var] <- factor(pre.rec[, var], 
                                        levels = levels(org.data))
          post.rec[, var] <- factor(post.rec[, var], 
                                         levels = levels(org.data))
        }
      } else if (is.integer(org.data[, var])) {
        
        q.engine[[var]][["discrete"]] <- discrete <- 1L
      } else { # logicals and small-domain numerics coerced to a factor
        
        org.data[, var] <- factor(org.data[, var])
        if (!is.null(pre.rec) && !is.null(post.rec)) {
          
          pre.rec[, var] <- factor(pre.rec[, var], 
                                        levels = levels(org.data))
          post.rec[, var] <- factor(post.rec[, var], 
                                        levels = levels(org.data))
        }
      }

      unique.values <- levels(org.data[, var])
      q.engine[[var]][["unique.values"]] <- unique.values

      org.data[, var] <- as.integer(org.data[, var]) +
                              runif(length(org.data[, var]), -0.5, 0.5)
      
      if (!is.null(pre.rec) && !is.null(post.rec)) {
        
        pre.rec[, var] <- as.integer(pre.rec[, var]) +
          runif(length(pre.rec[, var]), -0.5, 0.5)
        post.rec[, var] <- as.integer(post.rec[, var]) +
          runif(length(post.rec[, var]), -0.5, 0.5)
      }
    }

    # perform the Adaptation
    qr.data <- org.data[, c(var, curr.parents), drop = FALSE]

    object <- quant.method(qr.data, A.root = FALSE, base.ind, ...)
    q.engine[[var]][["object"]] <- object
    
    # ------------------------ #
    #  Recourse Data Learning  #
    # ------------------------ #
    
    if (!is.null(pre.rec) && !is.null(post.rec) && copula) {
      
      rec.idx <- is_updated(rec_act, var, adj.mat)
      
      if (sum(rec.idx) > 0) { # recourse data for this node exists
        
        empirical.pre <- computeQuants(
          q.engine[[var]][["object"]],
          pre.rec[rec.idx, c(var, curr.parents), drop = FALSE],
          pre.rec[rec.idx, curr.parents, drop = FALSE],
          rec.idx, test = TRUE, emp.only = TRUE, ...
        )
        
        val.pre <- pre.rec[, var] 
        U.pre <- vapply(seq_along(val.pre), 
                        function(x) ecdf(empirical.pre[x, ]) (val.pre[x]),
                        numeric(1L))
        
        empirical.post <- computeQuants(
          q.engine[[var]][["object"]],
          post.rec[rec.idx, c(var, curr.parents), drop = FALSE],
          post.rec[rec.idx, curr.parents, drop = FALSE],
          rec.idx, test = TRUE, emp.only = TRUE, ...
        )
        
        val.post <- post.rec[, var] 
        U.post <- vapply(seq_along(val.post), 
                         function(x) ecdf(empirical.post[x, ]) (val.post[x]),
                         numeric(1L))
        
        theta <- infer_theta(U.pre, U.post)
        if (theta == 35) {
          
          q.engine[[var]][["tau"]] <- 1
        } else if (theta == -35) {
          
          q.engine[[var]][["tau"]] <- -1
        } else {
          
          q.engine[[var]][["tau"]] <- BiCopPar2Tau(family = 5, par = theta)
        }
        
        # -------------------- #
        #  Hypothesis Testing  #
        # -------------------- #
        
        # are U.pre and U.post coming from a Frank copula?
        
        
        
        p.val <- frank_testing(U.pre, U.post, theta, 
                               browse = (var == "ExternalRiskEstimate"))
        q.engine[[var]][["p.val"]] <- p.val
        
        if (soft.test) {
          
          if (p.val < 0.01) message("Frank Copula rejected (soft test).\n")
        } else {
          
          assert_that(
            p.val > 0.01, 
            msg = paste0("Frank Copula hypothesis rejected for",
                         " variable ", var, ". \n",
                         "Supply `pre.rec` and `post.rec` arguments",
                         " without specifying `train.data`.")            
          )
        }
      }
      
    } else if (!is.null(pre.rec) && !is.null(post.rec) && !copula) {
      
      # -------------------- #
      #  No Copula Learning  #
      # -------------------- #
      
      post.dat <- post.rec[, c(var, curr.parents)]
      pre.dat <- pre.rec[, c(var, curr.parents)]
      names(post.dat) <- paste0(names(post.dat), "_post")
      names(pre.dat) <- paste0(names(pre.dat), "_pre")
      
      q.engine[[var]][["object.dual"]] <- ranger::ranger(
        as.formula(paste0(var, "_post ~ .")), data = cbind(pre.dat, post.dat),
        quantreg = TRUE, keep.inbag = TRUE, min.node.size = 20, ...
      )
      
      q.engine[[var]][["object.shadowdual"]] <- ranger::ranger(
        as.formula(paste0(var, "_post ~ .")), data = post.dat,
        quantreg = TRUE, keep.inbag = TRUE, min.node.size = 20, ...
      )
    }
  }

  structure(
    list(
      train = train.data,
      formula = formula,
      graph = ig,
      quant.method = deparse(substitute(quant.method)),
      adapt.call = match.call(),
      adj.mat = adj.mat,
      top.ord = top.ord,
      q.engine = q.engine
    ),
    class = "fairecourse"
  )
}

