\name{kullback}
\alias{kullback}
\title{Kullback-Leibler divergence}
\usage{
  kullback(proportion1, p1, mu1, proportion2, p2, mu2, m)
}
\arguments{
  \item{proportion1,proportion2}{a vector (which sum is
  equal 1) containing the proportion of the K clusters of
  the mixture}

  \item{p1,p2}{a matrix of size K*D, where K is the number
  of clusters and D the number of dimension, containing the
  probability of a good comparaison of the model}

  \item{mu1,mu2}{a matrix of size K*sum(m), containing the
  reference ranks. A row contains the reference rank for a
  cluster. In case of multivariate rank, for a cluster, the
  reference rank for each dimension are set successively on
  the same row.}

  \item{m}{a vector containing the size of rank for each
  dimension}
}
\value{
  a real, the Kullback-Leibler divergence
}
\description{
  This function compute the Kullback-Leibler divergence
  between 2 parametres of ISR for mixture of multivariate
  full rank
}
\examples{
proportion1=c(0.4,0.6)
p1=matrix(c(0.8,0.75),nrow=2)
mu1=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
proportion2=c(0.43,0.57)
p2=matrix(c(0.82,0.7),nrow=2)
mu2=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
dK=kullback(proportion1,p1,mu1,proportion2,p2,mu2,4)
}
\author{
  Quentin Grimonprez
}
\references{
  http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
}

