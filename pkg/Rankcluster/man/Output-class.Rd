\name{Output-class}
\alias{Output-class}
\title{Constructor of Output class}
\description{
  This class contains a result of a run.
}
\details{
  \describe{ \item{proportion}{a vector of size K (number
  of clusters), containing the estimates proportions of the
  mixture.} \item{pi}{a matrix of size K*p (p= number of
  dimensions), containing the estimates probabilities of
  good paired comparison of ISR (dispersion parameters).}
  \item{mu}{a matrix with K rows containing the estimates
  modal rankings of ISR.} \item{ll}{the estimated
  log-likelihood.} \item{bic}{the estimated BIC criterion.}
  \item{icl}{the estimated ICL criterion.} \item{tik}{a
  matrix of size n*K (n= number of individuals) containing
  the posterior probabilities for the ith observation to
  belong to the kth cluster.} \item{partition}{a vector of
  size n containing the estimated partition of the
  observations.} \item{entropy}{a matrix of size n*2
  containing the entropy of each individual and its
  estimated cluster.} \item{probability}{a matrix of size
  n*2 containing the probability of each individual and its
  estimated cluster.} \item{convergence}{if FALSE, no
  convergence has been achieved, no result is available.}
  \item{partial}{if FALSE, there is no partial rank in the
  data.} \item{partialRank}{a matrix containing the full
  ranking (estimated for the partial ranking and observed
  otherwise).} \item{distanceProp}{a list containing the
  distance between final proportion and proportion at each
  iteration of the algorithm.} \item{distancePi}{a list
  containing the distance between final probabilities of
  good paired comparison estimation and their current
  values at each iteration of the algorithm.}
  \item{distanceMu}{a list containing the distance between
  the final modal rankings estimation and their current
  values at each iteration of the algorithm.}
  \item{distanceZ}{a vector containing the Rand index
  between final partition and partition at each iteration
  of the algorithm.} \item{distancePartialRank}{a list
  containing the Kendall distance between final partial
  rankings estimation and their current values at each
  iteration of the algorithm.} \item{proportionInitial}{a
  vector containing the initialization of the proportions
  in the algorithm.} \item{piInitial}{a matrix containing
  the initialization of the probabilities of good paired
  comparison in the algorithm.} \item{muInitial}{a matrix
  containing the initialization of modal rankings in the
  algorithm.} \item{partialRankInitial}{a matrix containing
  the initialization of the partial rankings in the
  algorithm.} }
}

