
<!-- README.md is generated from README.Rmd. Please edit that file -->
VARCP
=====

``` r
# install.packages("devtools")
#devtools::install_github("MasimovR/VARCP")
```

``` r
library(VARCP)
#> Loading required package: MTS

phi <- matrix(c(0.2,-0.6,0.3,1),2,2)
sig <- matrix(c(4,0.8,0.8,1),2,2)
m1=VARMAsim(200, arlags = 1, phi=phi, sigma=sig)
m2=VARMAsim(200, arlags = 1, phi=phi, sigma=sig*2)
m <- rbind(m1$series, m2$series)  
model <- VAR(m, p = 1)
#> Constant term: 
#> Estimates:  -0.1353926 0.04900441 
#> Std.Error:  0.1231221 0.06312115 
#> AR coefficient matrix 
#> AR( 1 )-matrix 
#>        [,1]  [,2]
#> [1,]  0.220 0.321
#> [2,] -0.569 0.997
#> standard error 
#>        [,1]   [,2]
#> [1,] 0.0482 0.0431
#> [2,] 0.0247 0.0221
#>   
#> Residuals cov-mtx: 
#>          [,1]     [,2]
#> [1,] 5.754266 1.263857
#> [2,] 1.263857 1.512403
#>   
#> det(SSE) =  7.105433 
#> AIC =  1.98086 
#> BIC =  2.020774 
#> HQ  =  1.996666
```

``` r
CUSUM <- statcomp(model, type = "CUSUM")
str(CUSUM)
#> List of 2
#>  $ values: Named num [1:386] -0.206 -0.234 -0.275 -0.305 -0.339 ...
#>   ..- attr(*, "names")= chr [1:386] "8" "9" "10" "11" ...
#>  $ type  : chr "CUSUM"
#>  - attr(*, "class")= chr "VARCD"
```

``` r
plot(CUSUM, a = 0.01)
```

![](README-unnamed-chunk-3-1.png)

``` r
changetest(CUSUM)
#> 
#>  CUSUM test
#> 
#> data:  CUSUM
#> Estimated change point = 206, p-value = 3.127e-12
#> sample estimates:
#>     C(h) 
#> 3.686745
```

``` r
CUSUM2 <- statcomp(model, type = "CUSUM", trim = 50)
plot(CUSUM2)
```

![](README-unnamed-chunk-5-1.png)

``` r
LRT <- statcomp(model, type = "LRT")
DET <- statcomp(model, type = "DARLING-ERDOS")
```
