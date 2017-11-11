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
  class(statisitcvalue) <- "VARCD"
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
  class(statisitcvalue) <- "VARCD"
  return(statisitcvalue)
}


darling_erdos <- function(p, k, n, resi, S_t, points, X) {

  iSigma <- solve(S_t)

  S_h1 <- list()
  S_h2 <- list()
  sigmak <- list()
  for(i in points) {
    S_h1[[i]] <- MTS::VAR(X[1:i,, drop = FALSE], p = p,  output = F, include.mean = F)
    S_h2[[i]] <- MTS::VAR(X[(i + 1):n,, drop = FALSE], p = p,  output = F, include.mean = F)
  }

  prod3 <- function(X, order, Sigma, iSigma) {
    return( sum(diag(iSigma %*% Sigma * (nrow(X) - order) )))
  }

  QT <- 0
  QT2 <- 0
  a <- 0
  for(c in 1:nrow(resi)){
    a <- ((resi[c,,drop = F] %*% iSigma ) %*% resi[c, , drop = T])
    QT <- a + QT
    QT2 <- a^2 + QT2
  }

  Qk <- NULL
  QK <- NULL
  g <- 0
  for(c in points){
    g <- g + 1
    Qk[g] <- prod3(S_h1[[c]][["data"]], p, S_h1[[c]][["Sigma"]], iSigma)
    QK[g] <- prod3(S_h2[[c]][["data"]], p, S_h2[[c]][["Sigma"]], iSigma)
  }

  LR1 <- c(QT) - Qk - QK
  A <- NULL
  g <- 0
  for(i in points){
    g <- g + 1
    A[g] <- sum(diag(iSigma %*% S_h1[[i]][["Sigma"]] - diag(k)))^2
  }
  ks <- QT2/n
  LR2 <- LR1 + c( ( points - c(p)) /c(ks-k^2))*A
  d2 <- k*(k*p+1) + 1
  lgt <- 2*log(log(n))
  bt <- ((lgt + (d2/2)*log(log(log(n))) - lgamma(d2/2))^2) / lgt
  at <- sqrt(bt/lgt)
  values <- (LR2-bt)/at
  names(values) <- points
  statisitcvalue <- list("values" = values, type = "DARLING-ERDOS")
  class(statisitcvalue) <- "VARCD"
  return(statisitcvalue)
}

check <- function(x, model){
  if(is.null(model[[x]])) {
    stop(paste0("'",x,"'", " not found"))
  }
}
