# rbsrp
R package "rbsrp" provides a procedure to select response variables for ultrahigh dimensional multi-outcome data, where both the dimensions of response and predictor outcomes are substantially greater than the sample size. It also provides the multi-test procedure including Benjamini-Hochberg and Bonferroni correction.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/rbsrp")

# Usage

   - [x] [rbsrp-manual.pdf](https://github.com/xliusufe/rbsrp/inst/rbsrp-manual.pdf) ---------- Details of the usage of the package.

# Example
    library(rbsrp)
    library(stats)
    library(MASS)
    library(Matrix)

    n   <- 200
    p   <- 5
    q   <- 10
    q0  <- 5
    beta <- matrix(runif(p*q0),p,q0)
    eps <- matrix(rnorm(n*q),n,q)
    x <- matrix(rnorm(n*p),n,p)
    y <- cbind(x%*%beta, matrix(0,n,q-q0)) + eps

    fit <- rbsrp(y,x)
    fit$bestset

    fit <- bonf(y,x,alpha=0.05)
    fit$bestset

    fit <- bh(y,x,alpha=0.05)
    fit$bestset

    data(simulatedData)
    y = simulatedData$Y
    w = simulatedData$W
    fit <- rbsrp(y,w)
    fit$bestset
    
# References

Benjamini, Y. and Hochberg,  Y. (1995). Controlling the False Discovery Rate A Practical and Powerful Approach to Multiple testing. Journal of the Royal Statistical Society: Series B (Methodological). 57(1), 289-300.

Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. The Annals of Statistics, 29, 1165–1188.

Liu, X., Hu, J., and Liu, X. (2024). Random Projection-Based Response Best-Subset Selector for Ultrahigh Dimensional Multi-Outcome Data. Manuscript.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
