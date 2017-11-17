---
title: "VARCP"
output: html_document
---




```r
#devtools::install_github("MasimovR/VARCP")

library(VARCP)

phi <- matrix(c(0.2,-0.6,0.3,1),2,2)
sig <- matrix(c(4,0.8,0.8,1),2,2)
m1=VARMAsim(200, arlags = 1, phi=phi, sigma=sig)
m2=VARMAsim(200, arlags = 1, phi=phi, sigma=sig*2)
m <- rbind(m1$series, m2$series)  
model <- VAR(m, p = 1)
```

```
## Constant term: 
## Estimates:  0.1940187 0.09947162 
## Std.Error:  0.1159028 0.05881981 
## AR coefficient matrix 
## AR( 1 )-matrix 
##        [,1]  [,2]
## [1,]  0.123 0.244
## [2,] -0.617 0.970
## standard error 
##        [,1]   [,2]
## [1,] 0.0503 0.0422
## [2,] 0.0255 0.0214
##   
## Residuals cov-mtx: 
##          [,1]     [,2]
## [1,] 5.263437 1.057398
## [2,] 1.057398 1.355591
##   
## det(SSE) =  6.016978 
## AIC =  1.814585 
## BIC =  1.8545 
## HQ  =  1.830392
```


```r
CUSUM <- statcomp(model, type = "CUSUM")
str(CUSUM)
```

```
## List of 2
##  $ values: Named num [1:386] -0.242 -0.285 -0.256 -0.29 -0.324 ...
##   ..- attr(*, "names")= chr [1:386] "8" "9" "10" "11" ...
##  $ type  : chr "CUSUM"
##  - attr(*, "class")= chr "VARCD"
```

```r
plot(CUSUM, a = 0.01)
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-38-1.png)

```r
changetest(CUSUM)
```

```
## 
## 	CUSUM test
## 
## data:  CUSUM
## Estimated change point = 180, p-value = 7.394e-07
## sample estimates:
##     C(h) 
## 2.721268
```

```r
CUSUM2 <- statcomp(model, type = "CUSUM", trim = 50)
plot(CUSUM2)
```

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40-1.png)

```r
LRT <- statcomp(model, type = "LRT")
DET <- statcomp(model, type = "DARLING-ERDOS")
```
