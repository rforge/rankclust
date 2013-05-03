\name{Output-class}
\alias{Output-class}
\title{Constructor of Output class}
\description{
  This class contains a result of a run
}
\details{
  \describe{ \item{proportion}{A vector of size g (= the
  number of cluster), containing the estimates proportions
  of the mixture.} \item{pi}{a matrix of size g*d (d= the
  number of dimension), containing the estimates
  probabilities paramaeter of ISR.} \item{mu}{a matrix with
  g rows containing the estimates reference rank of ISR.}
  \item{ll}{The log-likelihood} \item{bic}{BIC criterion}
  \item{icl}{ICL criterion} \item{tik}{A matrix of size n*g
  (n=number of individuals) containing the posterior
  probabilities.} \item{partition}{A vector of size n
  containing for each indivuduals, the index of its
  cluster.} \item{entropy}{A matrix of size n*1 containing
  the entropy of each individuals.} \item{probability}{A
  matrix of size n*g containing the probability of each
  individuals in each cluster.} \item{convergence}{If
  FALSE, no convergence, no results available.}
  \item{partial}{If FALSE, there is no partial rank in the
  data.} \item{partialRank}{A matrix containing all the
  estimates of partial rank.} \item{distanceProp}{A matrix
  of size (Qsem-Bsem)*g containing the distance between
  final proportion and proportion at each iteration of the
  algorithm.} \item{distancePi}{A list containing the
  distance between final probabilities and probabilities at
  each iteration of the algorithm.} \item{distanceMu}{A
  list containing the distance between final mu and mu at
  each iteration of the algorithm.} \item{distanceZ}{A
  vector containing the Rand index between final partition
  and partition at each iteration of the algorithm.}
  \item{distancePartialRank}{A list containing the Kendall
  distance between final partial rank and partial rank at
  each iteration of the algorithm.}
  \item{proportionInitial}{A vector containing the
  initialization of the proportion in the algorithm.}
  \item{piInitial}{A matrix containing the initialization
  of tpi in the algorithm.} \item{muInitial}{A matrix
  containing the initialization of mu in the algorithm.}
  \item{partialRankInitial}{A matrix containing the
  initialization of the partial rank in the algorithm.} }
}

