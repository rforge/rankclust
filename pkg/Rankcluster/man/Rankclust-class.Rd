\name{Rankclust-class}
\alias{Rankclust-class}
\title{Constructor of Rankclust class}
\description{
  This class contains results of rankclust function.
}
\details{
  \describe{ \item{K}{list of character string with the
  estimation algorithm.  Possible values: "EM", "SEM",
  "CEM", c("EM","SEM"). Default value is "EM".}
  \item{data}{Data used in algorithm.}
  \item{criterion}{criterion defined in rankclust function
  to select the best result.} \item{convergence}{If 0, no
  convergence, no result available in results.}
  \item{results}{A list of the same length than K,
  containing Output objects.} }
}

