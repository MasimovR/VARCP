#' Plots estimated test statistic
#'
#' This function creates plot of computed test statistic with ability to add critical values as a red line
#'
#' @param VARCP A list of class "VARCP"
#' @param a Level of significance for which critical values are computed
#' @param ... Additional graphical arguments
#'
#' @return Displays a plot with automaticaly computed critical values
#'
#' @examples
#' ## Simulation of time series of length 200
#' phi <- matrix(c(0.2,-0.6,0.3,1),2,2)
#' sig <- matrix(c(4,0.8,0.8,1),2,2)
#' data <- VARMAsim(200, arlags = 1, phi = phi, sigma = sig)
#' model <- VAR(data[["series"]], p = 1)
#'
#' ## Test estimation
#' CUSUM <- statcomp(model, type = "CUSUM")
#'
#' Visualizing computed test statistic
#' plot(CUSUM, a = 0.1)
#' plot(CUSUM, a = 0.05)
#'
plot.VARCP <- function(VARCP, a = 0.05, ...) {
  if ((!is.numeric(a) | length(a) != 1)) stop('"a" must be numeric of length one')
  if (is.numeric(a) && (a < 0 | a >= 1)) stop('"a" must be positive and less than one')

  for (x in c("values", "type")) check(x, model = VARCP)


  x <- as.numeric(names(VARCP[["values"]]))
  y <- VARCP[["values"]]

  if(VARCP[["type"]] == "DET") {
    cvalue <- -2*log(-0.5*log(1-a))
    plot(ylim =c(min(y) , max(c(max(y), cvalue))), x = c(min(x), max(x)), y = rep(cvalue,2) , type = "l", col = "red",  main = "Darling-Erdös-type test statisic",
         xlab = "Time index", ylab = "Statistic value")
    lines(x, y)
  }

  if(VARCP[["type"]] == "CUSUM") {
    p <- 1 - a
    p2 <- 0
    cvalue <-0
    while(p != round(p2, 4)){
      cvalue <- cvalue + 0.0001
      p2 <- 1 + 2*sum( ((-1)^(1:250))*exp(-2*((1:250)^2)*cvalue^2))
    }
    plot(ylim =c(min(c(min(y), -cvalue)) , max(c(max(y), cvalue))), x = c(min(x), max(x)), y = rep(cvalue,2) , type = "l", col = "red",  main = paste0(VARCP[["type"]]," statistic"),
         xlab = "Time index", ylab = "Statistic value")
    lines(x = c(min(x), max(x)), y = rep(-cvalue,2) , type = "l", col = "red")
    lines(x, y)
  }
  if(VARCP[["type"]] == "LRT") {
    plot( x, y , type = "l",  main = paste0(VARCP[["type"]]," statistic"),
         xlab = "Time index", ylab = "Statistic value")
  }
}


#' Computes P value of estimated test statistic
#'
#' This function allows testing the presence of the change point by computing P value using corresponding asymptotic distribution of test statistic.
#'
#' @param VARCP A list of class "VARCP"
#'
#' @return changetest returns a list of class "htest" with following items
#' \item{method}{Character string which states type of the input test statistic}
#' \item{data.name}{Character string which states name of the input argument}
#' \item{p.value}{P value estimated using asymptotic distribution of the test statistic}
#' \item{estimate}{Maximal value of the test statistic}
#' \item{statistic}{Time index of estimated changepoint}
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
#' DET <- statcomp(model, type = "DET")
#'
#' ## Computing test P values
#' changetest(CUSUM)
#' changetest(DET)
#'
changetest <- function(VARCP) {

  for (i in c("values", "type")) check(i, VARCP)
  test <- list()
  test$method <- paste0(VARCP[["type"]], " test")
  test$data.name <- deparse(substitute(VARCP))

  if(VARCP[["type"]] == "CUSUM") {
    Cmax <- max(abs(VARCP[["values"]]))
    h <- as.numeric(names(which.max(abs(VARCP[["values"]])))) + 1
    attr(Cmax, "names") <- "C(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- 1 - (1 + 2*sum( ((-1)^(1:500))*exp(-2*((1:500)^2)*Cmax^2)))

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  if(VARCP[["type"]] == "DET") {
    Cmax <- max(VARCP[["values"]])
    h <- as.numeric(names(which.max(VARCP[["values"]]))) + 1
    attr(Cmax, "names") <- "G(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- 1 - (exp(-2*exp(-(Cmax/2))))

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  if(VARCP[["type"]] == "LRT") {
    Cmax <- max(abs(VARCP[["values"]]))
    h <- as.numeric(names(which.max(abs(VARCP[["values"]])))) + 1
    attr(Cmax, "names") <- "G(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- NULL

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  return(test)
}



