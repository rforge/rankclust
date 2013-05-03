# EMmeluni
# Fonction pour lancer EM unidimensionnel pour modeles de melange
#
EMmeluni<-function(X,g,maxIt=30,eps=1e-6,detail=FALSE)
{
	m=ncol(X)-1
	n=nrow(X)
	frequence=X[,m+1]
	rank=X[,-(m+1)]

	#Verification des donnees
	if(sum(apply(rank,1,checkRank,m))!=n)
		stop("Data are not correct")

	res=.Call("EMmelR",rank,frequence,m,g,eps,maxIt,detail,PACKAGE="Rankcluster")

	return(res)

}


# EMmelmulti
# fonction pour lancer EM multi dim pour modele de melange
EMmelmulti<-function(X,m,g,maxIt=30,eps=1e-6,detail=FALSE)
{
	n=nrow(X)
	d=length(m)
	if((ncol(X)-1)!=sum(m))
		stop("the number of column of data does not match to vector m")

	frequence=X[,sum(m)+1]
	rank=X[,-(sum(m)+1)]

	#Verification des donnees
	for(i in 1:d)
	{
		if(sum(apply(rank[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkRank,m[i]))!=n)
			stop("Data are not correct")
	}
	
	res=.Call("EMmelmultR",rank,frequence,m,g,eps,maxIt,detail,PACKAGE="Rankcluster")

	return(res)

}

# mixtureEM
# @title algorithm EM for 
# @author Grimonprez Quentin
# @param X matrix where each row is a rank and the last column contains the frequencies
# @param g number of groups
# @param m a vector with the size of ranks for each dimension
# @param maxIt the maximum number of iteration of the algorithm
# @param eps the threshold for loglikelihood convergency
# @note rank must be in ordering representation 
# @return a list containing the refererence rank mu, the probability pi of a correct comparison, the loglikelihood, the BIC and the ICL
# @references "Model-based clustering for rank data based on an insertion sorting algorithm", C. Biernacki, J. Jacques
# @examples
# data(quiz)#m=c(4,4,4,4)
# mixtureEM(quiz$frequency,K=2,m=quiz$m)
# @export

mixtureEM<-function(X,g,m=ncol(X)-1,maxIt=30,eps=1e-6,detail=FALSE)
{
	
	#algorithm
	if(length(m)==1)
	{
		if(m!=ncol(X)-1)
		{
			print(paste0("You put m=",m,", but X has ",ncol(X)," columns(rank of size ",ncol(X)-1," and 1 for the frequence)."))
 			print(paste0("The algorithm will continue with m=",ncol(X)-1))
		}

		res=EMmeluni(X,g,maxIt,eps,detail)
		res$referenceRank=liste2d2matG(res$referenceRank)
		res$pi=as.matrix(res$pi)
		colnames(res$pi)="dim 1"
		rownom=c()
		for(i in 1:g)	
			rownom=c(rownom,paste0("cl ",i))
		rownames(res$pi)=rownom
		
    
    res$partition=res$partition+1

    #entropie=[indice|entropie|classe|proba]
    index=c(1:nrow(X))
    #for(i in 1:nrow(X))
			#	index=c(index,rep(i,X[i,ncol(X)]))

    entropie=cbind(index,res$entropie,res$partition)
    maxProb=cbind(index,res$maxProb,res$partition)## a modif
    colnames(entropie)=c("index","entropy","cluster")		
    colnames(maxProb)=c("index","probability","cluster")


		result=new(Class="Output",bic=res$bic,icl=res$icl,
			LL=res$loglikelihood,proportion=res$proportion,pi=res$pi,
			mu=res$referenceRank,tik=res$tik,partition=res$partition,entropy=entropie,
			probability=maxProb,partial=FALSE)

		if(detail)
		{
			cat("RESULTS:\n")
			cat("NUMBER OF CLUSTERS: ",g)
			cat("\nLoglikelihood =",res$loglikelihood)   
			cat("\nBIC=",res$bic)
			cat("\nICL=",res$icl)		
			cat("\nProportion:",res$proportion)
			cat("\nProbabilities pi:\n")
			print(res$pi)
			cat("\nReference ranks mu:\n")
			print(res$referenceRank)
		}
	}
	else
	{
		if(length(m)>1)
		{
			res=EMmelmulti(X,m,g,maxIt,eps,detail)

			res$referenceRank=liste3d2mat(res$referenceRank)

			colnom=c()
			rownom=c()
			for(i in 1:g)	
				rownom=c(rownom,paste0("cl ",i))
			for(i in 1:length(m))	
				colnom=c(colnom,paste0("dim ",i))
			rownames(res$pi)=rownom
			colnames(res$pi)=colnom

			res$partition=res$partition+1
			#entropie=[indice|entropie|classe|proba]
			index=c(1:nrow(X))
			#for(i in 1:nrow(X))
			#	index=c(index,rep(i,X[i,ncol(X)]))
			entropie=cbind(index,res$entropie,res$partition)
			maxProb=cbind(index,res$maxProb,res$partition)
			colnames(entropie)=c("index","entropy","cluster")		
			colnames(maxProb)=c("index","probability","cluster")


			result=new(Class="Output",bic=res$bic,icl=res$icl,
				LL=res$loglikelihood,proportion=res$proportion,
				pi=res$pi,mu=res$referenceRank,tik=res$tik,
				partition=res$partition,entropy=entropie,
			probability=maxProb,partial=FALSE)

			if(detail)
			{
				cat("RESULTS:\n")
				cat("NUMBER OF CLUSTERS: ",g)
				cat("\nLoglikelihood =",res$loglikelihood)   
				cat("\nBIC=",res$bic)
				cat("\nICL=",res$icl)		
				cat("\nProportion:",res$proportion)
				cat("\nProbabilities pi:\n")
				print(res$pi)
				cat("\nReference ranks mu:\n")
				print(res$referenceRank)
			}
		}
		else
			stop("incorrect value of m")
	}

	return(result)
}

