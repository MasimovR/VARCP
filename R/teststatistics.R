CUSUM <- function(p, k, n, resi, S_t, points) {
  iSigma <- solve(S_t)
  Ah <- numeric(0)
  a <- 0
  for(c in 1:nrow(resi)){
    a <- ((resi[c,,drop = F] %*% iSigma ) %*% resi[c,]) + a
    Ah[c] <- a
  }
  A <- c(rep(NA, p), Ah)[points]
  values <- (points/sqrt(2*k*n)) * ((A/points) - c(a/n))
  names(values) <- points
  statisitcvalue <- list("values" = values, type = "CUSUM")
  class(statisitcvalue) <- "VARCP"
  return(statisitcvalue)
}

LRT <- function(p, k, n, resi, S_t, points) {

  ###
  v <- points/n
  vv <- 1-v
  ###
  S <- getDeterminant(S_t)
  S1 <- NA
  s1 <- 0
  for(t in 1:nrow(resi)){
    s1 <- crossprod(resi[t,, drop = F]) + s1
    S1[t] <- getDeterminant(s1/t)
  }
  S2 <- NA
  s2 <- 0
  for(t in nrow(resi):1){
    s2 <- crossprod(resi[t,, drop = F]) + s2
    S2[t] <- getDeterminant(s2/(n-t))
  }

  S1l <- (c(rep(NA, p), S1)[points])^v
  S2l  <- (c(rep(NA, p), S2)[points])^vv

  values <- n*log( S / (S1l * S2l) )
  names(values) <- points
  statisitcvalue <- list("values" = values, type = "LRT")
  class(statisitcvalue) <- "VARCP"
  return(statisitcvalue)
}


darling_erdos <- function(p, k, n, resi, S_t, points, X) {

  iSigma <- solve(S_t)

  X1 <- lapply(points, function (point) X[1:point,, drop = FALSE] )
  X2 <- lapply(points, function (point) X[(point + 1 - p):nrow(X),, drop = FALSE] )

  S_h1 <- lapply(X1,  MTS::VAR, p = p,  output = FALSE, include.mean = FALSE)
  S_h2 <- lapply(X2,  MTS::VAR, p = p,  output = FALSE, include.mean = FALSE)


  prod3 <- function(XN, iSigma) {
    Sigma <- XN[["Sigma"]]
    nx <- nrow(XN[["residuals"]])
    return( nx * sum(diag(iSigma %*% Sigma )))
  }

  QT <- 0
  QT2 <- 0
  a <- 0
  for(c in 1:nrow(resi)){
    a <- ((resi[c,,drop = F] %*% iSigma ) %*% resi[c,])
    QT <- a + QT
    QT2 <- a^2 + QT2
  }

  Qk <- sapply(S_h1, prod3, iSigma)
  QK <- sapply(S_h2, prod3, iSigma)

  LR1 <- c(QT) - Qk - QK

  A <- sapply(S_h1, function(sh) sum(diag( iSigma %*% sh[["Sigma"]] - diag(k) ))^2)

  ks <- QT2/2
  LR2 <- LR1 + c( ( points - c(p) ) / c(ks-k^2) )*A
  d2 <- k*(k*p+1) + 1
  lgt <- 2*log(log(n))
  bt <- ( (lgt + (d2/2)*log((log(log(n)))) - log(gamma(d2/2)))^2 ) / lgt
  at <- sqrt(bt/lgt)
  values <- (LR2-c(bt))/c(at)
  names(values) <- points
  statisitcvalue <- list("values" = values, type = "DET")
  class(statisitcvalue) <- "VARCP"
  return(statisitcvalue)
}

check <- function(x, model){
  if(is.null(model[[x]])) {
    stop(paste0("'",x,"'", " not found"))
  }
}
