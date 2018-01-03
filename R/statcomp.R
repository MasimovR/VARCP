#' Computes structural change test statistic
#'
#'
#' Computes test statistic value for each time point.
#' Likelihood-ratio test (LRT), Cumulative sum test (CUSUM) and Darling-Erdös-type test (DET) statistics are supported.
#'
#'
#' @param model A list returned by VAR or VARMA functions from MTS package
#' @param type A character string which determines type of the test statistic returned by the function, either "CUSUM", "LRT" or "DET". The default is "CUSUM"
#' @param trim A number of time points on data edges for which test statistic is not computed
#'
#' @return statcomp returns a list of class "VARCP" with following items
#' \item{values}{Named vector of values of selected test statistic, where name is time index for which statistic is computed}
#' \item{type}{Character string which states type of the computed test statistic}
#' @exportClass VARCD
#'
#'
#' @examples
#' ## Simulation of time series of length 200
#' phi <- matrix(c(0.2,-0.6,0.3,1),2,2)
#' sig <- matrix(c(4,0.8,0.8,1),2,2)
#' data <- VARMAsim(200, arlags = 1, phi = phi, sigma = sig)
#' model <- VAR(data[["series"]], p = 1)
#'
#' ## Tests estimation
#' CUSUM <- statcomp(model, type = "CUSUM")
#' LRT <- statcomp(model, type = "LRT", trim = 15)
#'

statcomp <- function(model, type = "CUSUM", trim = NULL) {
  is.list(model) || stop('"model" must be a list')

  if(is.null(model[["order"]]) & !is.null(model[["ARorder"]])) {
    model[["order"]] <- model[["ARorder"]]
  }

  for (x in c("data", "order", "residuals", "order", "Sigma")) check(x, model)

  if (!is.character(type) | length(type) != 1) stop('"type" must be character of length one')
  if (!(type %in% c("CUSUM", "LRT", "DET"))) stop('Invalid "type" argument specification')

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
  if(type == "DET") return(darling_erdos(p, k, n, resi, S_t, points, X))

}
