#####################################################################################################
# Response best-subst selector (RBS): # alpha in (0,0.05], gamma=0.95+2*alpha is from 0.95 to 1.05. 
#####################################################################################################
rbsrp<-function(Y,X,steplen=0.00125){
	n<-dim(Y)[1]
	q<-dim(Y)[2]
	p<-dim(X)[2]

	one<-matrix(1,n,1)
	projection1<-one%*%t(one)/n

	max_min = TRUE

	Y<-(diag(n)-projection1)%*%Y
	X<-(diag(n)-projection1)%*%X

	projection<-X%*%ginv(t(X)%*% X)%*% t(X)
	hat.sigma<-matrix(0,q,q)

	for (j in 1:q){
		hat.sigma[j,j]<-(t(Y[,j]) %*% (diag(n)-projection) %*% Y[,j]/(n-p))^(1/2) 
	}
	Y<-Y %*% solve(hat.sigma)

	bb<-floor(0.1/steplen)

	cv<-matrix(0,1,bb)
	gcv<-matrix(0,1,bb)
	alpha.set<-seq(0,0.1,steplen)

	delta<-matrix(0, bb, q)

	for (kk in 1:bb){    
		alpha<-alpha.set[kk+1]
		gamma<-0.95+2*alpha
		c<-qnorm(1-alpha, mean=0,sd=1,lower.tail = TRUE, log.p = FALSE)
		d<-p+c*sqrt(2*n*p/(n-p))     
		lambda<-(n-p)/d^(gamma)

		for (jj in 1:q) 
		{
			A<-t(Y[,jj]) %*% (diag(n)-projection) %*% Y[,jj]
			B<-(t(Y[,jj]) %*% projection %*% Y[,jj])^(gamma)
			if (A/B<=lambda) 
			{delta[kk,jj]<-1} 
		}

		Ps<-diag(c(delta[kk,]))
		residual<-(Y-projection%*%Y%*%Ps)

		J1<-c()
		for (j in 1:q)
		{
			if (delta[kk,j]==1){J1<-cbind(J1,j)}
		}

		dd<-n*q*(1-p*length(J1)/(n*q))^2

		CV.matrix<-matrix(0,n,q)
		for (ii in 1:n)
		{
			for (jj in 1:q)
			{
			CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-projection[ii,ii]*Ps[jj,jj])^2
			}
		}
		cv[kk]<-sum(CV.matrix)/(n*q)
		gcv[kk]<-norm((Y-projection%*%Y%*%Ps),"F")^2/dd
	}

	if (max_min==TRUE)
	{
		nnumber<-max(which.min(cv),which.min(gcv))
	}else{
		a1<-max(which.min(cv),which.min(gcv))+floor(0.01/steplen)
		if (a1<=bb){nnumber<-a1}else{nnumber<-bb}
	}

	decision.y<-delta[nnumber,]
	bestset.y<-c()

	for (j in 1:q){
		if (delta[nnumber,j]==1){bestset.y<-cbind(bestset.y,j)}
	} 

	result<-list(bestset=bestset.y, bestdecision=decision.y)
	return(result)
}