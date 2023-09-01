
correctInput <- function(formula, train.data, test.data, adj.mat, cfd.mat,
  top.ord, prot.attr, resolving.variables, quant.method) {

  if (is.null(adj.mat)) {

    assert_that(!is.null(top.ord))
    ap.nms <- top.ord

  } else {

    assert_that(
      matNames(adj.mat),
      withinRange(adj.mat),
      isAcyclic(adj.mat)
    )

    adj.mat <- adj.mat[colnames(adj.mat), ]
    ap.nms  <- colnames(adj.mat)
  }

  assert_that(
    all(resolving.variables %in% colnames(adj.mat)),
    prot.attr %in% ap.nms
  )

  train.data <- model.frame(formula, train.data)

  if (!is.null(test.data)) {
    test.data <- test.data[, colnames(train.data)[-1L]]
  }

  assert_that(sum(is.na(train.data)) == 0, nrow(train.data) > 0)

  if (!is.null(test.data)) {
    assert_that(sum(is.na(test.data)) == 0)
  }

  assert_that(
    all(colnames(train.data) %in% ap.nms),
    msg = "Train data has columns that do not appear in the adjacency matrix"
  )
  
  assert_that(
    length(unique(train.data[[prot.attr]])) == 2L,
    msg = "Protected attribute is not binary"
  )

  invisible(NULL)
}

withinRange <- function(mat) all(mat %in% c(0, 1))

isAcyclic <- function(mat) {

  num.walks <- mat

  for (i in seq_len(nrow(mat))) {
    num.walks <- mat + num.walks %*% mat
  }

  sum(diag(num.walks)) == 0
}

matNames <- function(mat) {

  assert_that(
    !is.null(colnames(mat)), !is.null(rownames(mat)),
    !anyDuplicated(colnames(mat)), !anyDuplicated(rownames(mat))
  )

  setequal(colnames(mat), rownames(mat))
}




