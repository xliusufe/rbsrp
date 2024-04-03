####################################################
# FDR, Benjamini-Hochberg method (JRBBS,1995)
####################################################

bh<-function(Y,X,alpha){
	n<-dim(X)[1]
	p<-dim(X)[2]
	q<-dim(Y)[2]
	one<-matrix(1,n,1)
	projection1<-one%*%t(one)/n

	X<-(diag(n)-projection1)%*%X
	Y<-(diag(n)-projection1)%*%Y
	PX<-X %*% ginv(t(X)%*% X) %*% t(X)

	J1<-c()
	TT<-matrix(0, 1, q)
	T2<-matrix(0, 1, q)
	pp<-matrix(0,1,q)
	for (j in 1:q) {
		TT[j]<-((n-p)/p)*t(Y[,j]) %*% PX %*% Y[,j]/(t(Y[,j]) %*% (diag(n)-PX) %*% Y[,j])-1
		T2[j]<-TT[j]/sqrt(2*n/(p*(n-p)))
		pp[j]<-pnorm(T2[j], mean=0,sd=1,lower.tail =FALSE, log.p = FALSE)
	}
	OrderP<-sort(pp)
	Orderindex<-order(pp)
	K0<-c()
	for (i in 1:q) {
		if (OrderP[i]<= i*alpha/q) {
			K0<-c(K0, i)
		} else {
			K0<-K0
		}
	}
	if (length(K0)>0) {
		J1<-Orderindex[1:K0[length(K0)]];
		J1<-sort(J1)
	} else {
		J1<-J1
	} 
	result<-list(bestset=J1)
	return(result)
}


