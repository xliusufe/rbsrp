# First install the following packages:
# install.packages("devtools")
# library(devtools)
# install_github("xliusufe/rbsrp")
# install.packages("stats")
# install.packages("MASS")
# install.packages("Matrix")
library(rbsrp)
library(stats)
library(MASS)
library(Matrix)

assig <- function(n_args) {
    cargs <- vector("list", length(n_args))
    for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
    t(expand.grid(cargs))
}

gram.orth <- function(U){
    pp  <- dim(U)[2]
    V   <- matrix(0, pp, pp)
    PC  <- matrix(0, pp, pp)
    for (s in 1:pp)
    {
        W       <- PC %*% U[,s]
        V[,s]   <- (U[,s]-W)/norm((U[,s]-W),"2")
        if (s < pp)
        {
            PC <- PC + V[,s] %*% t(V[,s])
        }
    }
    return(V)
}

# Generate data
dataGenerator <- function(X, C, n, q, Sigma05, Coefficient, normal){
    if (normal == TRUE){
        Error0 <- matrix(rnorm(n*q,0,1),n,q) 
    }
    else{
        Error0 <- matrix(rt(n*q,3,ncp=0),n,q)
    }
    
    Error <- Error0 %*% Sigma05
    
    Y <- X %*% Coefficient + Error
    
    W <- X %*% C

    return(
        list(
            Y = Y,
            W = W
        )

    )
}

main <- function(n,p,q,m0){

    if(m0 == 100 && n == 100){theta.h = 2.4}
    else if(m0 == 100 && n == 200){theta.h = 1.5}
    else if(m0 == 1 && n == 100 && p == 1000){theta.h = 10}
    else if(m0 == 1 && n == 100 && p == 3000){theta.h = 19}
    else if(m0 == 1 && n == 200 && p == 1000){theta.h = 4.65}
    else if(m0 == 1 && n == 200 && p == 3000){theta.h = 8}

    nreplication    = 500                           # replication times
    q1              = 100
    normal          = TRUE                          # normal=TRUE means N(0,1)  and normal=FALSE means student t(3)

    Z  <- matrix(rnorm(n*p, 0, 1), n, p)
	p0 <- p/m0
	Q  <- matrix(0, p, p)
	for (i in 1:m0){
		U           <-matrix(rnorm(p0*p0, 0, 2), p0, p0)
		startpoint  <-(i-1)*p0 + 1
		endpoint    <-i*p0
		Q[startpoint:endpoint,startpoint:endpoint] <- gram.orth(U)
	}
	s  <- ceiling(n^(0.8))
	d1 <- matrix(1, 1, s)
	d2 <- matrix(0, 1, p-s)
	a0 <- 0
	for (i in 1:p-s){
		a0 <- a0 + 1/i^4
	}
	for (i in 1:p-s){
		d2[i] <- sqrt((n-s)*i^(-4)/a0)
	}
	dd <- c(d1, d2)
	D  <- diag(dd)

	X <- Z %*% Q %*% D %*% t(Q)

    p1 <- ceiling(n/2)
	x.activeset <- c(1:p1)

	k0 <- ceiling(0.49*n)
	C <- matrix(rnorm(p*k0, 0, 1), p, k0)

    y.activeset <- c(1:q1)
    y.inactiveset <- setdiff(c(1:q), y.activeset)

    Coefficient <- matrix(0, p, q)
    Theta0 <- matrix(runif(p1*q1, min=-theta.h, max=theta.h), p1, q1)
    Coefficient[x.activeset, y.activeset] <- Theta0

    MS1     <- matrix(0, 1, nreplication)
    TP1     <- matrix(0, 1, nreplication)
    FP1     <- matrix(0, 1, nreplication)
    Prec1   <- matrix(0, 1, nreplication)
    Fscore1 <- matrix(0, 1, nreplication)

    MS2     <-matrix(0, 1, nreplication)
    TP2     <-matrix(0, 1, nreplication)
    FP2     <-matrix(0, 1, nreplication)
    Prec2   <-matrix(0, 1, nreplication)
    Fscore2 <-matrix(0, 1, nreplication)

    MS3     <-matrix(0, 1, nreplication)
    TP3     <-matrix(0, 1, nreplication)
    FP3     <-matrix(0, 1, nreplication)
    Prec3   <-matrix(0, 1, nreplication)
    Fscore3 <-matrix(0, 1, nreplication)

    MS4     <-matrix(0, 1, nreplication)
    TP4     <-matrix(0, 1, nreplication)
    FP4     <-matrix(0, 1, nreplication)
    Prec4   <-matrix(0, 1, nreplication)
    Fscore4 <-matrix(0, 1, nreplication)

    sigma   <- sample(c(1:15), q, replace=TRUE)
    U1      <- matrix(rnorm(q*q, 0, 2), q, q)
    O1      <- gram.orth(U1)
    Sigma05 <- O1 %*% diag(sigma) %*% t(O1)

    ptm <- proc.time()
    # Loop of replications    
    for (m in 1:nreplication){
        dgp = dataGenerator(X, C, n, q, Sigma05, Coefficient, normal)
        Y   = dgp$Y
        W   = dgp$W
        
        ########################################################### 
        #  Method 1: RBS via max(min(cv),min(gcv))
        ###########################################################

        RBS         <- rbsrp(Y, W)
        rbs.bestset <- RBS$bestset
        
        MS1[m]  <- length(rbs.bestset)
        TP1[m]  <- length(intersect(y.activeset, rbs.bestset))
        FP1[m]  <- length(intersect(y.inactiveset, rbs.bestset))
        if (MS1[m]!=0){
            Prec1[m]    <- TP1[m]/MS1[m]
            Fscore1[m]  <- (2/q1)*Prec1[m]*TP1[m]/(Prec1[m]+TP1[m]/q1)
        }
        else{
            Prec1[m]    <-0
            Fscore1[m]  <-0
        }
        
        # ##########################################################################
        #  Method 2: Bonferoni
        # ##########################################################################
        
        Bonf.test    <- bonf(Y, W, 0.05)
        Bonf.bestset <- Bonf.test$bestset

        MS2[m]  <- length(Bonf.bestset)
        TP2[m]  <- length(intersect(y.activeset, Bonf.bestset))
        FP2[m]  <- length(intersect(y.inactiveset, Bonf.bestset))
        if (MS2[m]!=0){
            Prec2[m]    <- TP2[m]/MS2[m]
            Fscore2[m]  <- (2/q1)*Prec2[m]*TP2[m]/(Prec2[m]+TP2[m]/q1)
        }
        else{
            Prec2[m]    <- 0
            Fscore2[m]  <- 0
        }
        
        # ########################################################################
        # #  Method 3: FDR method
        # ########################################################################
        
        FDR05           <- bh(Y, W, 0.05)
        fdr05.bestset   <- FDR05$bestset
        
        MS3[m]  <-length(fdr05.bestset)
        TP3[m]  <-length(intersect(y.activeset, fdr05.bestset))
        FP3[m]  <-length(intersect(y.inactiveset, fdr05.bestset))
        if (MS3[m]!=0)
        {
            Prec3[m]    <- TP3[m]/MS3[m]
            Fscore3[m]  <- (2/q1)*Prec3[m]*TP3[m]/(Prec3[m]+TP3[m]/q1)
        }
        else{
            Prec3[m]    <-0
            Fscore3[m]  <-0
        }  
        
        # ########################################################################
        # #  Method 4: BY method
        # ########################################################################
        
        alpha1  <- 0.05/sum(c(1/q:1))
        
        BY          <- bh(Y, W, alpha1)
        by.bestset  <- BY$bestset

        MS4[m]  <- length(by.bestset)
        TP4[m]  <- length(intersect(y.activeset, by.bestset))
        FP4[m]  <- length(intersect(y.inactiveset, by.bestset))
        if (MS4[m]!=0)
        {
            Prec4[m]    <- TP4[m]/MS4[m]
            Fscore4[m]  <- (2/q1)*Prec4[m]*TP4[m]/(Prec4[m]+TP4[m]/q1)
        }
        else{
            Prec4[m]    <- 0
            Fscore4[m]  <- 0
        }
    }

    # Means of replication
    TPR1        <- mean(TP1, na.rm = TRUE)/q1
    FPR1        <- mean(FP1, na.rm = TRUE)/(q-q1)
    F1          <- mean(Fscore1, na.rm = TRUE)
    P1          <- mean(Prec1, na.rm = TRUE)
    RBS.mean    <- c(TPR1, 1-FPR1, F1)
    
    TPR2        <- mean(TP2, na.rm = TRUE)/q1
    FPR2        <- mean(FP2, na.rm = TRUE)/(q-q1)
    F2          <- mean(Fscore2, na.rm = TRUE)
    P2          <- mean(Prec2, na.rm = TRUE)
    Bonf.mean   <- c(TPR2, 1-FPR2, F2)
    
    TPR3        <- mean(TP3, na.rm = TRUE)/q1
    FPR3        <- mean(FP3, na.rm = TRUE)/(q-q1)
    F3          <- mean(Fscore3, na.rm = TRUE)
    P3          <- mean(Prec3, na.rm = TRUE)
    FDR05.mean  <- c(TPR3, 1-FPR3, F3)
    
    TPR4        <- mean(TP4, na.rm = TRUE)/q1
    FPR4        <- mean(FP4, na.rm = TRUE)/(q-q1)
    F4          <- mean(Fscore4, na.rm = TRUE)
    P4          <- mean(Prec4, na.rm = TRUE)
    BY.mean     <- c(TPR4, 1-FPR4, F4)

    result      = matrix(NA, 4, 3)
    result[1,]  = round(RBS.mean, digits=4)
    result[2,]  = round(Bonf.mean, digits=4)
    result[3,]  = round(FDR05.mean, digits=4)
    result[4,]  = round(BY.mean, digits=4)
    
    # ################

    print(result)
    cat("cpu time=", proc.time()-ptm, "\n")
    return(result)
}

exam <- function(){

    sample.size         = c(100,200)          # Sample size
    response.number     = c(2000,5000)        # The number of response variables
    predictor.number    = c(1000,3000)        # The number of predictor variables
    m00                 = c(100,1)            # two cases: m0=1 and m0=100

    n_args  = c(length(response.number), length(sample.size), length(m00), length(predictor.number)) # 2*2*2*2=16
    jobs    = assig(n_args)

    result = matrix(NA, 32, 6)
 
    for(number in 1:ncol(jobs)){
        id      = jobs[, number]           # number = 1: q=2000,n=100,m0=100,p=1000; number = 2: q=5000,n=100,m0=100,p=1000; number = 3: q=2000,n=200,m0=100,p=1000; number = 4: q=5000,n=200,m0=100,p=1000.
        q       = response.number[id[1]]
        n       = sample.size[id[2]]
        m0      = m00[id[3]]
        p       = predictor.number[id[4]]

        cat("n = ", n, " p = ", p, " q = ", q, " m0 = ", m0) 
        if(1 <= number && number <= (ncol(jobs)/2)){result[(4*number-3):(4*number),1:3] = main(n,p,q,m0) }
        else {result[(4*number-35):(4*number-32),4:6] = main(n,p,q,m0)}
    }
    colnames(result) = rep(c("TPR", "TNR", "F-score"), 2)
    rownames(result) = rep(c("RBS", "Bonf", "FDR05", "BY"), 8)
    write.csv(result, file = "./Table1.csv")
}

exam()