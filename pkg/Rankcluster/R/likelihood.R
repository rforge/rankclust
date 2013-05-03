#' This function estimates the loglikelihood, bic and icl.
#' @title criteria estimation
#' @author Quentin Grimonprez
#' @param data a matrix which each row is a rank (partial or not; for partial rank, missing elements of a rank are put to 0 ) (the last column contains frequencies of each rank).
#' @param proportion a vector (which sum is equal 1) containing the proportion of the K clusters of the mixture
#' @param p a matrix of size K*D, where K is the number of clusters and D the number of dimension, containing the probability of a good comparaison of the model
#' @param mu a matrix of size K*sum(m), containing the reference ranks. A row contains the reference rank for a cluster. In case of multivariate rank, for a cluster, the reference rank for each dimension are set successively on the same row.
#' @param m a vector containing the size of rank for each dimension
#' @param Ql number of iterations of the Gibbs sampler for estimation of log-likelihood (only for SEM algorithm, default value=100)
#' @param Bl burn-in period for estimation of log-likelihood (only for SEM algorithm, default value=50)
#' @return a list containing
#' \item{ll}{loglikelihood estimation}
#' \item{bic}{bic criterion}
#' \item{icl}{icl criterion}
#' @examples
#' data(quiz)
#' res=rankclust(quiz$data,m=quiz$m,K=2)
#' crit=criteria(quiz$data,res[2]@@proportion,res[2]@@pi,res[2]@@mu,quiz$m)
#' @export
criteria <-function(data,proportion,p,mu,m,Ql=500,Bl=100)
{
  if(missing(proportion))
    stop("proportion is missing")
  if(missing(mu))
    stop("mu is missing")
  if(missing(p))
    stop("p is missing")
  if(missing(m))
    stop("m is missing")

  
  #data
  if(missing(data))
    stop("data is missing")
  if(!is.numeric(data) || !is.matrix(data))
    stop("data must be a matrix of positive integer")
  if(length(data[data>=0])!=length(data))
    stop("data must be a matrix of positive integer")

  #proportion
  if(!is.vector(proportion,mode="numeric"))
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(min(proportion)<0)
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(abs(1-sum(proportion))>1e-10)
    stop("proportion must be a vector of positive real whose sum equal 1")
  
  
  #m
  if(!is.vector(m,mode="numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(length(m)!=length(m[m>1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(!min(m==round(m)))
    stop("m must be a (vector of) integer strictly greater than 1")
  if( (length(m)!=ncol(p)) )
    stop("The number of column of p and m don't match.")
  if(sum(m)!=ncol(mu)) 
    stop("The number of column of mu and sum(m) don't match.")

  #p
  if(!is.numeric(p) || !is.matrix(p))
    stop("p must be a matrix of probabilities")
  if( (min(p)<0) && (max(p)>1) )
    stop("p must be a matrix of probabilities")
  if( (nrow(p)!=length(proportion)) || (nrow(p)!=nrow(mu)) )
    stop("The number of rows of p doesn't match with the others parameters.")
  
  #Ql
  if(!is.numeric(Ql) || (length(Ql)>1))
    stop("Ql must be a strictly positive integer")
  if( (Ql!=round(Ql)) || (Ql<=0))
    stop("Ql must be a strictly positive integer")
  
  #Bl
  if(!is.numeric(Bl) || (length(Bl)>1))
    stop("Bl must be a strictly positive integer lower than Ql")
  if( (Bl!=round(Bl)) || (Bl<=0) || (Bl>=Ql))
    stop("Bl must be a strictly positive integer lower than Ql")
  
  for(i in 1:length(m))
  {
    if(sum(apply(data[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkPartialRank,m[i]))!=nrow(data))
      stop("Data are not correct")
  }

  
  a=t(p)
  
  LL=.Call("loglikelihood",data,mu,a,proportion,m,Ql,Bl,PACKAGE="Rankcluster")
  
  return(LL)
}


