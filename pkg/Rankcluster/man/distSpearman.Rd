\name{distSpearman}
\alias{distSpearman}
\title{Spearman's distance between two ranks}
\usage{
  distSpearman(x, y)
}
\arguments{
  \item{x,y}{two ranks of size m}
}
\value{
  an integer, the Spearman distance between x and y
}
\description{
  Spearman's distance is a square of Euclidean distance
  between two rank vector
}
\examples{
x=1:5
y=c(2,3,1,4,5)
distSpearman(x,y)
}
\author{
  Julien Jacques
}

