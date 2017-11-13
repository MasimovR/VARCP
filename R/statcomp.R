statcomp <- function(model = list(...), type = "CUSUM", trim = NULL) {
  is.list(model) || stop('"model" must be a list')

  if(is.null(model[["order"]]) & !is.null(model[["ARorder"]])) {
    model[["order"]] <- model[["ARorder"]]
  }

  for (x in c("data", "order", "residuals", "order", "Sigma")) check(x, model)

  if (!is.character(type) | length(type) != 1) stop('"type" must be character of length one')
  if (!(type %in% c("CUSUM", "LRT", "DARLING-ERDOS"))) stop('Invalid "type" argument specification')

  if ((!is.numeric(trim) | length(type)  != 1) & !is.null(trim)) stop('"trim" must be numeric of length one')
  if (is.numeric(trim) && trim < 0) stop('"trim" must be positive')

  X <- model[["data"]]
  p <- model[["order"]]
  iSigma <- solve(model[["Sigma"]])

  resi <- model[["residuals"]]
  S_t <- model[["Sigma"]]

  k <- ncol(X)
  n <- nrow(X)
  if (is.null(trim)) trim <- k*(p + 1) + k + 1
  points <- (trim + 1):(n-trim)

  if(type == "CUSUM") return(CUSUM(p, k, n, resi, S_t, points))
  if(type == "LRT") return(LRT(p, k, n, resi, S_t, points))
  if(type == "DARLING-ERDOS") return(darling_erdos(p, k, n, resi, S_t, points, X))

}
