plot.VARCD <- function(VARCD, a = 0.05) {
  if ((!is.numeric(a) | length(a) != 1)) stop('"a" must be numeric of length one')
  if (is.numeric(a) && (a < 0 | a >= 1)) stop('"a" must be positive and less than one')

  for (x in c("values", "type")) check(x, model = VARCD)


  x <- as.numeric(names(VARCD[["values"]]))
  y <- VARCD[["values"]]

  if(VARCD[["type"]] == "DARLING-ERDOS") {
    cvalue <- -2*log(-0.5*log(1-a))
    plot(ylim =c(min(y) , max(c(max(y), cvalue))), x = c(min(x), max(x)), y = rep(cvalue,2) , type = "l", col = "red",  main = "Darling-Erdös-type test statisic",
         xlab = "Time index", ylab = "Statistic value")
    lines(x, y)
  }

  if(VARCD[["type"]] == "CUSUM") {
    p <- 1 - a
    p2 <- 0
    cvalue <-0
    while(p != round(p2, 4)){
      cvalue <- cvalue + 0.0001
      p2 <- 1 + 2*sum( ((-1)^(1:250))*exp(-2*((1:250)^2)*cvalue^2))
    }
    plot(ylim =c(min(c(min(y), -cvalue)) , max(c(max(y), cvalue))), x = c(min(x), max(x)), y = rep(cvalue,2) , type = "l", col = "red",  main = paste0(VARCD[["type"]]," statistic"),
         xlab = "Time index", ylab = "Statistic value")
    lines(x = c(min(x), max(x)), y = rep(-cvalue,2) , type = "l", col = "red")
    lines(x, y)
  }
  if(VARCD[["type"]] == "LRT") {
    plot( x, y , type = "l",  main = paste0(VARCD[["type"]]," statistic"),
         xlab = "Time index", ylab = "Statistic value")
  }
}

changetest <- function(x) {

  for (i in c("values", "type")) check(i, x)
  test <- list()
  test$method <- paste0(x[["type"]], " test")
  test$data.name <- deparse(substitute(x))

  if(x[["type"]] == "CUSUM") {
    Cmax <- max(abs(x[["values"]]))
    h <- as.numeric(names(which.max(abs(x[["values"]])))) + 1
    attr(Cmax, "names") <- "C(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- 1 - (1 + 2*sum( ((-1)^(1:250))*exp(-2*((1:250)^2)*Cmax^2)))

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  if(x[["type"]] == "DARLING-ERDOS") {
    Cmax <- max(abs(x[["values"]]))
    h <- as.numeric(names(which.max(abs(x[["values"]])))) + 1
    attr(Cmax, "names") <- "G(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- 1 - (exp(-2*exp(-(Cmax/2))))

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  if(x[["type"]] == "LRT") {
    Cmax <- max(abs(x[["values"]]))
    h <- as.numeric(names(which.max(abs(x[["values"]])))) + 1
    attr(Cmax, "names") <- "G(h)"
    attr(h, "names") <- "Estimated change point"
    test$p.value <- NULL

    test$estimate <- Cmax
    test$statistic <- h

    class(test) <- "htest"
  }
  return(test)
}



