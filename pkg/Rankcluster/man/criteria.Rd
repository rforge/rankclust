\name{criteria}
\alias{criteria}
\title{criteria estimation}
\usage{
  criteria(data, proportion, p, mu, m, Ql = 500, Bl = 100)
}
\arguments{
  \item{data}{a matrix which each row is a rank (partial or
  not; for partial rank, missing elements of a rank are put
  to 0 ) (the last column contains frequencies of each
  rank).}

  \item{proportion}{a vector (which sum is equal 1)
  containing the proportion of the K clusters of the
  mixture}

  \item{p}{a matrix of size K*D, where K is the number of
  clusters and D the number of dimension, containing the
  probability of a good comparaison of the model}

  \item{mu}{a matrix of size K*sum(m), containing the
  reference ranks. A row contains the reference rank for a
  cluster. In case of multivariate rank, for a cluster, the
  reference rank for each dimension are set successively on
  the same row.}

  \item{m}{a vector containing the size of rank for each
  dimension}

  \item{Ql}{number of iterations of the Gibbs sampler for
  estimation of log-likelihood (only for SEM algorithm,
  default value=100)}

  \item{Bl}{burn-in period for estimation of log-likelihood
  (only for SEM algorithm, default value=50)}
}
\value{
  a list containing \item{ll}{loglikelihood estimation}
  \item{bic}{bic criterion} \item{icl}{icl criterion}
}
\description{
  This function estimates the loglikelihood, bic and icl.
}
\examples{
data(quiz)
res=rankclust(quiz$data,m=quiz$m,K=2)
crit=criteria(quiz$data,res[2]@proportion,res[2]@pi,res[2]@mu,quiz$m)
}
\author{
  Quentin Grimonprez
}

