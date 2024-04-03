# First install the following packages:
# install.packages("devtools")
# library(devtools)
# install_github("xliusufe/rbsrp")
# install.packages("stats")
# install.packages("MASS")
# install.packages("Matrix")
library(rbsrp)
library(stats)
library(Matrix)
library(MASS)

rbsrp_1 <- function(X, Y, nmethods, ntest, p, q, alpha){
	indy_rbs 	= rep(0, q)
	if(nmethods == 1){#  Method 1: RBS via max(min(cv),min(gcv))
        RBS     <- rbsrp(Y, X)
		indy 	<- RBS$bestset
	}
	else if(nmethods == 2){#  Method 2: Bonferoni
        Bonf.test   <- bonf(Y, X, 0.05)
		indy 		<- Bonf.test$bestset
	}
	else if(nmethods == 3){#  Method 3: FDR method
        FDR05   <- bh(Y, X, 0.05)
		indy 	<- FDR05$bestset
	}
	else if(nmethods == 4){#  Method 4: BY method
		alpha1	<-0.05/sum(c(1/q:1))
        BY      <- bh(Y, X, alpha1)
		indy	<- BY$bestset
	}
	indy_rbs[indy] = 1
	return( indy_rbs )
}

main_rbss <- function(nchrom, NS, methodj, len_sul = NS, id_sul = 1){
	alpha   <- 0.05
	Sigma2  <- NULL
	Tn 		<- NULL
	RSS     <- 0
	PE 		<- 0
	Est_rbs <- NULL

	# load real data ###
	load("breastdata.rda")
	X0	= t(breastdata$dna[breastdata$chrom==nchrom,])          # CNV
	Y	= t(breastdata$rna[which(breastdata$genechr==nchrom),]) # gene

	q = dim(Y)[2]
	p = dim(X0)[2]
	n = dim(Y)[1]
	Y = scale(Y)

	k0	<- ceiling(0.49*n)
	if(p<k0){
		X = X0
	}
	else{
		if(methodj==1){set.seed(1)}
		C 	<- matrix(rnorm(p*k0,0,1),p,k0)
		X	<- X0%*%C
	}

	id_beg = (id_sul-1)*len_sul+1
	id_end = id_sul*len_sul
	for(j in id_beg:id_end){
		print(j)
		seed_id = 1e4+j
		set.seed(seed_id)
		cv.id    = sample.int(n,size = 10)
		Ytrain   = Y[-cv.id,]
		Xtrain   = X[-cv.id,]
		Ytest    = Y[cv.id,]
		Xtest    = X[cv.id,]

		indy     <- rbsrp_1(Xtrain, Ytrain, methodj, n-length(cv.id), p, q, alpha)
		Est_rbs  <- cbind(Est_rbs, indy)

		NPx = solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)
		IPx = diag(n-length(cv.id)) - Xtrain%*%NPx

		nonzero.loc = which(indy==1)
		if(length(nonzero.loc)>0){
			RSS = RSS + sum((IPx%*%Ytrain[,nonzero.loc])^2)
			PE  = PE  + sum((Ytest[,nonzero.loc] - Xtest%*%NPx%*%Ytrain[,nonzero.loc])^2)
		}

	}
	return(list(RSS = RSS/len_sul/10/q, PE = PE/len_sul/10/q, Est_rbs = Est_rbs))
}

nchrom   = 3
NS       = 100
ptm 	 <- Sys.time()
ntop 	 = 80

RSS     	<- NULL
PE     		<- NULL
ngene10 	<- NULL
genename80 	= NULL
genepos10  	= NULL
nselect    	= NULL
ngenes_all 	= NULL
Pos_NT     	= NULL
result = matrix(NA, 80, 12)

for(methodj in 1:4){
	fit_rbss <- main_rbss(nchrom, NS, methodj, NS, 1)

	activeA = fit_rbss$Est_rbs
	fit_rbss$RSS
	fit_rbss$PE

	RSS = c(RSS, fit_rbss$RSS)
	PE  = c(PE,  fit_rbss$PE)

	geneid_sort = sort(rowSums(activeA),decreasing = T, index.return = T)
	ngene10 	= c(ngene10, sum(rowSums(activeA)>10))
	ngenes_all 	= c(ngenes_all, mean(colSums(activeA)))

	nselect = rbind(nselect, rowSums(activeA)[geneid_sort$ix[1:ntop]])

	load("breastdata.rda")
	genenames21 = breastdata$genenames[which(breastdata$genechr==nchrom)]
	genenames80 = genenames21[geneid_sort$ix[1:ntop]]
	id_na 		= which(genenames80=="")
	genenames80[id_na] = "--"
	genepos21   = breastdata$genepos[which(breastdata$genechr==nchrom)]
	Poss_NT = rbind(genepos21[geneid_sort$ix[1:ntop]], rowSums(activeA)[geneid_sort$ix[1:ntop]])
	result[,(methodj*3-2):(methodj*3)]=cbind((unname(genenames80)),t(Poss_NT))
}
cat("nchrom", nchrom, '\n')
print(rbind(RSS,PE, ngenes_all, ngene10))

cat("Sys time = ", format(Sys.time()-ptm),"\n")

colnames(result) = c("Gene name_rbs", "Position_rbs", "Times_rbs", "Gene name_bonf", "Position_bonf", "Times_bonf", "Gene name_bh", "Position_bh", "Times_bh", "Gene name_by", "Position_by", "Times_by")
write.csv(result, file = paste("./Table4.csv"))