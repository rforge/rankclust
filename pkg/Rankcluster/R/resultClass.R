###################################################################################
##' Constructor of Output class
##'
##' This class contains a result of a run
##'
##' \describe{
##'   \item{proportion}{A vector of size g (= the number of cluster), containing the estimates proportions of the mixture.}
##'   \item{pi}{a matrix of size g*d (d= the number of dimension), containing the estimates probabilities paramaeter of ISR.}
##'   \item{mu}{a matrix with g rows containing the estimates reference rank of ISR.}
##'   \item{ll}{The log-likelihood}
##'   \item{bic}{BIC criterion}
##'   \item{icl}{ICL criterion}
##'   \item{tik}{A matrix of size n*g (n=number of individuals) containing the posterior probabilities.}
##'   \item{partition}{A vector of size n containing for each indivuduals, the index of its cluster.}
##'   \item{entropy}{A matrix of size n*1  containing the entropy of each individuals.}
##'   \item{probability}{A matrix of size n*g containing the probability of each individuals in each cluster.}
##'   \item{convergence}{If FALSE, no convergence, no results available.}
##'   \item{partial}{If FALSE, there is no partial rank in the data.}
##'   \item{partialRank}{A matrix containing all the estimates of partial rank.}
##'   \item{distanceProp}{A matrix of size (Qsem-Bsem)*g containing the distance between final proportion and proportion at each iteration of the algorithm.}
##'   \item{distancePi}{A list containing the distance between final probabilities and probabilities at each iteration of the algorithm.}
##'   \item{distanceMu}{A list containing the distance between final mu and mu at each iteration of the algorithm.}
##'   \item{distanceZ}{A vector containing the Rand index between final partition and partition at each iteration of the algorithm.}
##'   \item{distancePartialRank}{A list containing the Kendall distance between final partial rank and partial rank at each iteration of the algorithm.}
##'   \item{proportionInitial}{A vector containing the initialization of the proportion in the algorithm.}
##'   \item{piInitial}{A matrix containing the initialization of tpi in the algorithm.}
##'   \item{muInitial}{A matrix containing the initialization of mu in the algorithm.}
##'   \item{partialRankInitial}{A matrix containing the initialization of the partial rank in the algorithm.}
##' }
##'
##'
##' @name Output-class
##' @rdname Output-class
## @exportClass Output
##'
setClass(
  Class="Output",
  representation=representation(
    proportion="numeric",
    pi="matrix",
    mu="matrix",
    ll="numeric",
    bic="numeric",
    icl="numeric",
    tik="matrix",
    partition="numeric",
	entropy="matrix",
	probability="matrix",
	convergence="logical",
	partial="logical",
	partialRank="matrix",
	distanceProp="list",
	distancePi="list",
	distanceMu="list",
	distanceZ="numeric",
	distancePartialRank="list",
	proportionInitial="numeric",
    piInitial="matrix",
    muInitial="matrix",
	partialRankInitial="matrix"
    ),
  prototype=prototype(
    proportion=numeric(0),
    pi=matrix(nrow=0,ncol=0),
    mu=matrix(nrow=0,ncol=0),
    ll=numeric(0),
    bic=numeric(0),
    icl=numeric(0),
    tik=matrix(nrow=0,ncol=0),
    partition=numeric(0),
	entropy=matrix(nrow=0,ncol=0),
	probability=matrix(nrow=0,ncol=0),
	convergence=logical(0),
	partial=logical(0),
	partialRank=matrix(nrow=0,ncol=0),
	distanceProp=list(),
	distancePi=list(),
	distanceMu=list(),
	distanceZ=numeric(0),
	distancePartialRank=list(),
	proportionInitial=numeric(0),
    piInitial=matrix(nrow=0,ncol=0),
    muInitial=matrix(nrow=0,ncol=0),
	partialRankInitial=matrix(nrow=0,ncol=0)
  )
)



###################################################################################
##' Constructor of Rankclust class
##'
##' This class contains results of rankclust function.
##'
##' \describe{
##'   \item{K}{list of character string with the estimation algorithm.  Possible values: "EM", "SEM", "CEM", c("EM","SEM"). Default value is "EM".}
##'   \item{data}{Data used in algorithm}
##'   \item{criterion}{criterion defined in rankclust function to select the best result.}
##'   \item{convergence}{If 0, no convergence, no results available in results.}
##'   \item{results}{A list of the same length than K, containing Output objects.}
##' }
##'
##'
##' @name Rankclust-class
##' @rdname Rankclust-class
## @exportClass Rankclust
##'
setClass(
  Class="Rankclust",
  representation=representation(
    K="numeric",
    results="list",
	data="matrix",
    criterion="character",
	convergence="logical"
  ),
	prototype=prototype(
    results=list(),
	data=matrix(ncol=0,nrow=0),
    K=numeric(0),
    criterion="bic",
	convergence=logical(0)    
  )

)

setMethod(
  f="[",
  signature="Rankclust",
  definition=function(x,i,j,drop){
	if(x@convergence)
	{
		if(is.numeric(i))
		{
			if(i %in% x@K)
			{
				return(x@results[[which(x@K==i)]])       
			}
			else
				stop("Invalid number of cluster.") 
		}
		else
		{
			if(i=="bic")
			{
				bic=rep(NA,length(x@K))
				for(iter in 1:length(x@K))
				{
					if(x@results[[iter]]@convergence)
						bic[iter]=x@results[[iter]]@bic
				}
				return(bic)
			}
			else
			{
				if(i=="icl")
				{
					icl=rep(NA,length(x@K))
					for(iter in 1:length(x@K))
					{
						if(x@results[[iter]]@convergence)
							icl[iter]=x@results[[iter]]@icl
					}
					return(icl)
				}
				else
				{
					if(i=="ll")
					{
						ll=rep(NA,length(x@K))
						for(iter in 1:length(x@K))
						{
							if(x@results[[iter]]@convergence)
								ll[iter]=x@results[[iter]]@ll
						}
						return(ll)
					}
					else
					{
						stop("Invalid Name.")
					}
		
				}

			}
		}
    }
  }
)

setMethod(
  f="summary",
  signature = "Rankclust",
  definition = function(object,...) {
	if(object@convergence)
	{
		if (object@criterion=="bic") 
		{
			BIC=c()
			for(i in object@K)
			{
				BIC=c(BIC,object@results[[which(object@K==i)]]@bic)
			}
			index=which(BIC==min(BIC))
		}
		else
		{
			ICL=c()
			for(i in object@K)
			{
				ICL=c(ICL,object@results[[which(object@K==i)]]@icl)
			}
			index=which(ICL==min(ICL))

		}

		cat("******************************************************************\n")
		cat("NUMBER OF CLUSTERS: ",object@K[index],"\n")
		if(object@criterion=="bic")  
			cat(object@criterion,"=",object[object@K[index]]@bic)
		else
			cat(object@criterion,"=",object[object@K[index]]@icl)
		cat("\nLoglikelihood =",object[object@K[index]]@ll)
		cat("\n\n*************************PARAMETERS*******************************\n")
		cat("Proportion:",object[object@K[index]]@proportion)
		cat("\nProbabilities pi:\n")
		print(object[object@K[index]]@pi)
		cat("\nReference ranks mu:\n")
		print(object[object@K[index]]@mu)
		cat("\n*************************CLUSTERING*******************************\n")
		cat("Ranks with the highest entropy for each cluster:\n")
		for(i in 1:object@K[index])
		{
			classe=object[object@K[index]]@entropy[object[object@K[index]]@entropy[,3]==i,]
#print(dim(classe))
#print(nrow(classe))
#print(nrow(object@data))##pb si rien et voir entropie
			if(length(classe)!=0)
			{
				if(length(classe)==3)
				{
					best5=classe
					print(cbind(object@data[best5[1],-ncol(object@data)],best5[2:3]))
				}
				else
				{
					best5=classe[order(classe[,2],decreasing=TRUE),][1:min(5,nrow(classe)),]
					print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
				}
			}
	
		}
		rm(best5)	
		cat("Ranks with the highest probability for each cluster:\n")
		for(i in 1:object@K[index])
		{
			classe=object[object@K[index]]@probability[object[object@K[index]]@probability[,3]==i,]
			if(length(classe)!=0)
			{
				if(length(classe)==3)
				{
					best5=classe
					print(cbind(object@data[best5[1],-ncol(object@data)],best5[2:3]))
				}
				else
				{
					best5=classe[order(classe[,2],decreasing=TRUE),][1:min(5,nrow(classe)),]
					print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
				}
				
			}	
		}   
	  	rm(classe)
		rm(best5)
		if(object[object@K[index]]@partial)
		{
			cat("\n*************************PARTIAL RANK*****************************\n")
			print(object[object@K[index]]@partialRank)
		}
		
		cat("\n******************************************************************\n")
	}
	else
		cat("\nNo convergence. Please retry\n")

  }  
)

setMethod(
  f="show",
  signature = "Output",
  definition = function(object) {
    cat("ll=",object@ll)
    cat("\nbic =",object@bic)
    cat("\nicl =",object@icl)
    cat("\nproportion:",object@proportion)
    cat("\npi:\n")
    print(object@pi)
    cat("\nmu:\n")
    print(object@mu)
    cat("\npartition:\n")
    print(object@partition[1:min(50,length(object@partition))])
    if(min(50,length(object@partition))==50)
      cat("\nOnly the first 50 are printed, total length:",length(object@partition))
    cat("\ntik:\n")
    print(object@tik[1:min(50,nrow(object@tik)),])
    if(min(50,nrow(object@tik))==50)
      cat("\nOnly the first 50 rows are printed, total rows:",nrow(object@tik))
    }
)


setMethod(
  f="show",
  signature = "Rankclust",
  definition = function(object) {
    for(i in object@K)
	{
		cat("\n******************************************************************\n")
		cat("Number of clusters:",i)
		cat("\n******************************************************************\n")
		show(object@results[[which(object@K==i)]])
		cat("\n******************************************************************\n")
		
	}
    }
)
