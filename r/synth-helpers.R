
reg_fun_beta <- function(var, top.ord, data) {
  
  idx <- which(top.ord == var)
  form_str <- paste(
    var, "~", paste(top.ord[seq_len(idx - 1)], collapse = "+"),
    "|", paste(top.ord[seq_len(idx - 1)], collapse = "+")
  )
  
  betareg(as.formula(form_str), data = data)
}

reg_fun_geom <- function(var, top.ord, data) {
  
  idx <- which(top.ord == var)
  form_str <- paste(
    var, "~", paste(top.ord[seq_len(idx - 1)], collapse = "+")
  )

  glm(as.formula(form_str), data, family = negative.binomial(theta = 1))
}

reg_fun <- function(var, top.ord, data, type, ...) {
  
  if (type == "beta") return(reg_fun_beta(var, top.ord, data))
  if (type == "geom") return(reg_fun_geom(var, top.ord, data))
}

gen_fun <- function(x, data, quants, ...) {
  
  UseMethod("gen_fun", x)
}

gen_fun.betareg <- function(x, data, quants) {
  
  alphbet <- beta_mom2par(predict(x, newdata = data, type = "response"), 
                          predict(x, newdata = data, type = "variance"))
  
  qbeta(quants, alphbet$alpha, alphbet$beta)
}

gen_fun.glm <- function(x, data, quants) {
  
  # obtain the distribution parameters
  p <- geom_mom2par(predict(x, newdata = data, type = "response"))
  
  qgeom(quants, p)
}


