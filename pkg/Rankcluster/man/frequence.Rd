\name{frequence}
\alias{frequence}
\title{Convert data}
\usage{
  frequence(X, m = ncol(X))
}
\arguments{
  \item{X}{a matrix containing ranks}

  \item{m}{a vector with the size of rank of each
  dimension}
}
\value{
  a matrix with rank and frequencies
}
\description{
  This function take in input a matrix containing all the
  ranks (a rank can be repeated) and returns a matrix which
  the m first column are the differents ranks observed and
  the last column contains the frequency of each rank.
}
\examples{
X=matrix(1:4,ncol=4,nrow=5,byrow=TRUE)
Y=frequence(X)
Y
}
\author{
  Quentin Grimonprez
}

