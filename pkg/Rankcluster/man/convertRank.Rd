\name{convertRank}
\alias{convertRank}
\title{change the representation of a rank}
\usage{
  convertRank(x)
}
\arguments{
  \item{x}{a rank of size m}
}
\value{
  a vector with the rank in the other representation
}
\description{
  convertRank convert a rank from ranking to ordering
  representation or ordering to ranking representation. The
  transformation to convert a rank fron ordering to ranking
  representation is the same that from ranking to ordering
  representation, there is no need to precise the
  representation of rank x.
}
\details{
  A rank in ordering representation contains first the
  number of the object ranks first, second the number of
  the object ranks second,....

  A rank in ranking representation contains first the rank
  assigned to the first object, second the rank assigned to
  the second object,....

  Let us consider the following example to illustrate both
  notations: a judge, which has torank three holidays
  destinations according to its preferences, O1 =
  Countryside, O2 =Mountain and O3 = Sea, ranks first Sea,
  second Countryside, and last Mountain. The ordering
  result of the judge is x = (3, 1, 2) whereas the ranking
  result is (2, 3, 1).
}
\examples{
x=c(2,3,1,4,5)
convertRank(x)
}
\author{
  Julien Jacques
}

