#'convertRank convert a rank from ranking to ordering representation or
#'  ordering to ranking representation.
#'The transformation to convert a rank fron ordering to ranking representation is the same that from ranking to ordering representation, there is no need to precise the representation of rank x.
#'
#'
#'A rank in ordering representation contains first the number of the object ranks first, second the number of the object ranks second,....
#'
#'
#'A rank in ranking representation contains first the rank assigned to the first object, second the rank assigned to the second object,....
#'
#'
#'Let us consider the following example to illustrate both notations: a judge, which has torank three holidays destinations according to its preferences, O1 = Countryside, O2 =Mountain and O3 = Sea, ranks first Sea, second Countryside, and last Mountain. The ordering result of the judge is x = (3, 1, 2) whereas the ranking result is (2, 3, 1).
#' @useDynLib Rankcluster
#' @title change the representation of a rank
#' @author Julien Jacques
#' @param x a rank of size m 
#' @return a vector with the rank in the other representation
#' @examples
#' x=c(2,3,1,4,5)
#' convertRank(x)
#' @export

convertRank <- function(x){
	return(sort.int(x,index.return=1)$ix)	
}

# checkRank  check if a vector is a rank


checkRank <- function(x,m=length(x))
{
	if(sum(sort(x)==(1:m))==m)
		return(TRUE)
	else
		return(FALSE)	
}

# checkPartialRank check if a vector is a rank


checkPartialRank <- function(x,m=length(x))
{
	if((length(x[x<=m])==m)&& (length(x[x>=0])==m) && (length(unique(x[x!=0]))==length(x[x!=0])))
		return(TRUE)
	else
		return(FALSE)
}

# completeRank complete partial that have only one missing element


completeRank <-function(x)
{
	if(length(x[x==0])==1)
	{	
		m=length(x)
		a=1:m
		a[x[x!=0]]=0
		x[x==0]=a[a!=0]
	}
	return(x)
}

#' This function take in input a matrix containing all the ranks (a rank can be repeated) and returns a matrix which the m first column are the differents ranks observed and the last column contains the frequency of each rank.
#' @title Convert data
#' @author Quentin Grimonprez
#' @param X a matrix containing ranks
#' @param m a vector with the size of rank of each dimension
#' @return a matrix with rank and frequencies
#' @examples
#' X=matrix(1:4,ncol=4,nrow=5,byrow=TRUE)
#' Y=frequence(X)
#' Y
#' @export
frequence <-function(X,m=ncol(X))
{
	if(missing(X))
		stop("X is missing")
	if(!is.numeric(X) || !is.matrix(X))
		stop("X must be a matrix of positive integer")
	if(length(X[X>=0])!=length(X))
		stop("X must be a matrix of positive integer")
	if(!is.vector(m,mode="numeric"))
		stop("m must be a (vector of) integer strictly greater than 1")
	if(length(m)!=length(m[m>1]))
		stop("m must be a (vector of) integer strictly greater than 1")

	if(length(m)==1)
	{
		if(m!=ncol(X))
		{
			print(paste0("You put m=",m,", but X has ",ncol(X)," columns(rank of size ",ncol(X)-1," and 1 for the frequence)."))
 			print(paste0("The algorithm will continue with m=",ncol(X)-1))
		}
	}

	res=.Call("freqMultiR",X,m,PACKAGE="Rankcluster")
	
	data=matrix(0,ncol=length(res$data[[1]])+1,nrow=length(res$data))
	for(i in 1:nrow(data))
		data[i,]=c(res$data[[i]],res$freq[[i]])
	

	return(data)

}

#' This function simulate a sample of univariate ranks according to the ISR(p,mu)
#' @title simulate a sample of ISR(p,mu)
#' @author Julien Jacques
#' @param n size of the sample
#' @param p probability of correct comparaison according to mu
#' @param mu reference rank in ordering representation
#' @return a matrix with simulated ranks
#' @references 
#' [1] C.Biernacki and J.Jacques (2012), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, in press, DOI 10.1016/j.csda.2012.08.008
#' @examples
#' x=simulISR(30,0.8,1:4)
#' @export
simulISR <-function(n,p,mu)
{
	if(missing(n))
		stop("n is missing")
	if(missing(mu))
		stop("mu is missing")
	if(missing(p))
		stop("p is missing")

	if(!is.numeric(n) || (length(n)>1))
		stop("n must be a strictly positive integer")
	if( (n!=round(n)) || (n<=0))
		stop("n must be a strictly positive integer")

	if(!is.numeric(p) || (length(n)>1))
		stop("p must be a real between 0 and 1")
	if( (p>1) || (p<0))
		stop("p must be a real between 0 and 1")

	if(!is.vector(mu,mode="numeric"))
		stop("mu must be a complete rank")
	if(!checkRank(mu))
		stop("mu must be a complete rank")



	res=.Call("simulISRR",n,length(mu),mu,p,PACKAGE="Rankcluster")

	return(res)
}

