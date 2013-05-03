#' This function compute the Kullback-Leibler divergence between 2 parametres of ISR for mixture of multivariate full rank
#' @title Kullback-Leibler divergence 
#' @author Quentin Grimonprez
#' @param proportion1,proportion2 a vector (which sum is equal 1) containing the proportion of the K clusters of the mixture
#' @param p1,p2 a matrix of size K*D, where K is the number of clusters and D the number of dimension, containing the probability of a good comparaison of the model
#' @param mu1,mu2 a matrix of size K*sum(m), containing the reference ranks. A row contains the reference rank for a cluster. In case of multivariate rank, for a cluster, the reference rank for each dimension are set successively on the same row.
#' @param m a vector containing the size of rank for each dimension
#' @return a real, the Kullback-Leibler divergence 
#' @references 
#' http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
#' @examples
#' proportion1=c(0.4,0.6)
#' p1=matrix(c(0.8,0.75),nrow=2)
#' mu1=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' proportion2=c(0.43,0.57)
#' p2=matrix(c(0.82,0.7),nrow=2)
#' mu2=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' dK=kullback(proportion1,p1,mu1,proportion2,p2,mu2,4)
#' @export
kullback <-function(proportion1,p1,mu1,proportion2,p2,mu2,m)
{
	if(missing(proportion1))
		stop("proportion1 is missing")
	if(missing(mu1))
		stop("mu1 is missing")
	if(missing(p1))
		stop("p1 is missing")
	if(missing(proportion2))
		stop("proportion2 is missing")
	if(missing(mu2))
		stop("mu2 is missing")
	if(missing(p2))
		stop("p2 is missing")

	eps=1e-10
	#proportion1
	if(!is.vector(proportion1,mode="numeric"))
		stop("proportion1 must be a vector of positive real whose sum equal 1")
	if(min(proportion1)<0)
		stop("proportion1 must be a vector of positive real whose sum equal 1")
	if(abs(1-sum(proportion1))>eps)
		stop("proportion1 must be a vector of positive real whose sum equal 1")

	#proportion2
	if(!is.vector(proportion2,mode="numeric"))
		stop("proportion2 must be a vector of positive real whose sum equal 1")
	if(min(proportion2)<0)
		stop("proportion2 must be a vector of positive real whose sum equal 1")
	if(abs(1-sum(proportion2))>eps)
		stop("proportion2 must be a vector of positive real whose sum equal 1")
  if(length(proportion1)!=length(proportion2))
    stop("proportion1 and proportion2 must have the same length.")
  
	#m
	if(!is.vector(m,mode="numeric"))
	  stop("m must be a (vector of) integer strictly greater than 1")
	if(length(m)!=length(m[m>1]))
	  stop("m must be a (vector of) integer strictly greater than 1")
	if(!min(m==round(m)))
	  stop("m must be a (vector of) integer strictly greater than 1")
	if( (length(m)!=ncol(p1)) || (length(m)!=ncol(p2)) )
	  stop("The number of column of p1 or p2 and m don't match.")
	if( (sum(m)!=ncol(mu1)) || (sum(m)!=ncol(mu2)) )
	  stop("The number of column of mu1 or mu2 and sum(m) don't match.")
  
	#p1
	if(!is.numeric(p1) || !is.matrix(p1))
	  stop("p1 must be a matrix of probabilities")
	if( (min(p1)<0) && (max(p1)>1) )
	  stop("p1 must be a matrix of probabilities")
	
	#p2
	if(!is.numeric(p2) || !is.matrix(p2))
	  stop("p2 must be a matrix of probabilities")
	if( (min(p2)<0) && (max(p2)>1) )
	  stop("p2 must be a matrix of probabilities")
	if(length(p1)!=length(p2))
	  stop("p1 and p2 must have the same size.")
  if( (nrow(p1)!=length(proportion1)) || (nrow(p1)!=nrow(mu1)) )
    stop("The number of rows of p1 doesn't match with the others parameters.")
  if( (nrow(p2)!=length(proportion2)) || (nrow(p2)!=nrow(mu2)) )
    stop("The number of rows of p2 doesn't match with the others parameters.") 

	a=t(p1)
	b=t(p2)
	dKL=.Call("kullback",m,mu1,mu2,a,b,proportion1,proportion2,PACKAGE="Rankcluster")

	return(dKL)
}



#' This function compute the p-value of the khi2 adequation test (only for univariate data)
#' @title khi2 test
#' @author Quentin Grimonprez
#' @param data a matrix where each row contains a rank of size m
#' @param proportion a vector (which sum is equal 1) containing the proportion of the K clusters of the mixture
#' @param p a vector of size K, where K is the number of clusters, containing the probability of a good comparaison of the model
#' @param mu a matrix of size K*m, where m is the lengh of a rank, containing the reference ranks for the model
#' @param nBoot number of iteration of the bootstrap for estimating the p-value
#' @return a real, the p-value of the khi2 adequation test 
#' @examples
#' proportion=c(0.4,0.6)
#' p=c(0.8,0.75)
#' mu=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' data=rbind(simulISR(proportion[1]*100,p[1],mu[1,]),simulISR(proportion[2]*100,p[2],mu[2,]))
#' pval=khi2(data,proportion,mu,p)
#' @export
khi2 <-function(data,proportion,mu,p,nBoot=1000)
{
	if(missing(proportion))
		stop("proportion is missing")
	if(missing(mu))
		stop("mu is missing")
	if(missing(p))
		stop("p is missing")
	if(missing(data))
		stop("data is missing")

	eps=1e-10
	#proportion
	if(!is.vector(proportion,mode="numeric"))
		stop("proportion must be a vector of positive real whose sum equal 1")
	if(min(proportion)<0)
		stop("proportion must be a vector of positive real whose sum equal 1")
	if(abs(1-sum(proportion))>eps)
		stop("proportion must be a vector of positive real whose sum equal 1")

	#p
	if(!is.vector(p,mode="numeric"))
	  stop("p must be a vector of probabilities")
	if( (min(p)<0) && (max(p)>1) )
	  stop("p must be a vector of probabilities")
  
	#mu
	if(!is.numeric(mu) || !is.matrix(mu))
	  stop("mu must be a matrix of positive integer")
	if(min(mu)<1)
	  stop("mu must be a matrix of positive integer")
  if(nrow(mu)!=length(proportion))
    stop("The number of rows of mu and the length of proportion don't match.")
	if(nrow(mu)!=length(p))
	  stop("The number of rows of mu and the length of p don't match.")

	#data
	if(missing(data))
	  stop("data is missing")
	if(!is.numeric(data) || !is.matrix(data))
	  stop("X must be a matrix of positive integer")
	if(length(data[data>=0])!=length(data))
	  stop("data must be a matrix of positive integer")
  
  if(ncol(data)!=ncol(mu))
    stop("mu and data must have the same number of columns.")
	

	#nBoot
	if(!is.numeric(nBoot) )
    stop("nBoot must be a positive integer.")
	if(length(nBoot)!=1 )
	  stop("nBoot must be a positive integer.")
	if( (nBoot<0) || (nBoot!=round(nBoot)) )
	  stop("nBoot must be a positive integer.")


	pval=.Call("adkhi2partial",data,p,proportion,mu,nBoot,PACKAGE="Rankcluster")

	return(pval)
}
