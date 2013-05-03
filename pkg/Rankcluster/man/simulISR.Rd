\name{simulISR}
\alias{simulISR}
\title{simulate a sample of ISR(p,mu)}
\usage{
  simulISR(n, p, mu)
}
\arguments{
  \item{n}{size of the sample}

  \item{p}{probability of correct comparaison according to
  mu}

  \item{mu}{reference rank in ordering representation}
}
\value{
  a matrix with simulated ranks
}
\description{
  This function simulate a sample of univariate ranks
  according to the ISR(p,mu)
}
\examples{
x=simulISR(30,0.8,1:4)
}
\author{
  Julien Jacques
}
\references{
  [1] C.Biernacki and J.Jacques (2012), A generative model
  for rank data based on sorting algorithm, Computational
  Statistics and Data Analysis, in press, DOI
  10.1016/j.csda.2012.08.008
}

