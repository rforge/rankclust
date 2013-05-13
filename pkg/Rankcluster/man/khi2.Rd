\name{khi2}
\alias{khi2}
\title{khi2 test}
\usage{
  khi2(data, proportion, mu, pi, nBoot = 1000)
}
\arguments{
  \item{data}{a matrix where each row contains a rank of
  size m}

  \item{proportion}{a vector (which sum is equal 1)
  containing the proportion of the K clusters of the
  mixture}

  \item{pi}{a vector of size K, where K is the number of
  clusters, containing the probability of a good paired
  comparaison of the model (dispersion parameter)}

  \item{mu}{a matrix of size K*m, where m is the size of a
  rank, containing the modal rankings of the model
  (position parameter)}

  \item{nBoot}{number of bootstrap iterations used to
  estimate the khi2 adequation test p-value}
}
\value{
  a real, the p-value of the khi2 adequation test
}
\description{
  This function compute the p-value of the khi2 adequation
  test (only for univariate data)
}
\examples{
proportion=c(0.4,0.6)
pi=c(0.8,0.75)
mu=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
data=rbind(simulISR(proportion[1]*100,pi[1],mu[1,]),simulISR(proportion[2]*100,pi[2],mu[2,]))#simulate a data set with declared parameters.
pval=khi2(data,proportion,mu,pi)
}
\author{
  Quentin Grimonprez
}

