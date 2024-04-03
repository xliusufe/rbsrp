###################################################
# Bonferroni Test
###################################################
bonf<-function(Y, X, alpha){
	n<-dim(Y)[1]
	q<-dim(Y)[2]
	p<-dim(X)[2]

	one<-matrix(1,n,1)
	projection1<-one%*%t(one)/n

	X<-(diag(n)-projection1)%*%X
	Y<-(diag(n)-projection1)%*%Y

	PX<-X %*%ginv(t(X)%*% X)%*% t(X)
	c<-qnorm(1-alpha/q, mean=0,sd=1,lower.tail = TRUE, log.p = FALSE)
	d<-p+c*sqrt(2*n*p/(n-p))     
	lambda<-(n-p)/d

	delta<-matrix(0, 1, q)
	J1<-c()
	for (j in 1:q) {
		A<-t(Y[,j]) %*% (diag(n)-PX) %*% Y[,j]
		B<-(t(Y[,j]) %*% PX %*% Y[,j])
		if (A/B<=lambda) {delta[j]<-1; J1<-c(J1, j)} else {J1<-J1}
	}

	result<-list(decision.set=delta, bestset=J1)
	return(result)
}