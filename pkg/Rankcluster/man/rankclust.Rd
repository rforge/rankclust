\name{rankclust}
\alias{rankclust}
\title{model-based clustering for multivariate partial ranking}
\usage{
  rankclust(data, m = ncol(data), K = 1, criterion = "bic",
    algorithm = "SEM", Qsem = 100, Bsem = 20,
    RjSE = m * (m - 1)/2, RjM = m * (m - 1)/2, Ql = 500,
    Bl = 100, maxItEM = 30, epsEM = 1e-06, maxTry = 3,
    run = 1, detail = FALSE)
}
\arguments{
  \item{data}{a matrix which each row is a rank (partial or
  not; for partial rank, missing elements of a rank are put
  to 0 ) (the last column contains frequencies of each
  rank).}

  \item{m}{a vector with the size of the rank of each
  dimension (default value is the number of column of the
  data minus 1 (corresponding to the frequency)}

  \item{K}{an integer or a vector of integer with the
  number of group for clustering}

  \item{criterion}{"bic" or "icl", criterion to minize for
  selecting the number of clusters}

  \item{algorithm}{type of algorithm : "EM" (only for
  complete rank) or "SEM" for SEM-Gibbs algorithm(default
  value)}

  \item{Qsem}{the total number of iterations for the SEM
  algorithm (defaut value=40)}

  \item{Bsem}{burn-in period for SEM algorithm (default
  value=10)}

  \item{RjSE}{a vector containing the number of iteration
  for each dimension of the Gibbs algorithm in the SE step
  for generate partial ranks and orders of
  presentation(only for SEM algorithm, default
  value=m(m-1)/2)}

  \item{RjM}{a vector containing the number of iterations
  for each dimension for the Gibbs Sampler in the M
  step(only for SEM algorithm, default value=m(m-1)/2)}

  \item{Ql}{number of iterations of the Gibbs sampler for
  estimation of log-likelihood (only for SEM algorithm,
  default value=100)}

  \item{Bl}{burn-in period for estimation of log-likelihood
  (only for SEM algorithm, default value=50)}

  \item{maxItEM}{the maximum number of iterations of the EM
  algorithm (defaut value=30)}

  \item{epsEM}{the threshold for loglikelihood convergency
  (only for EM algorithm, defaut value=1e-6)}

  \item{maxTry}{number of try for the convergency of the
  SEM-Gibbs algorithm (default value=3)}

  \item{run}{number of runs of the algorithm for each value
  of K}

  \item{detail}{boolean, if TRUE, time and others
  informations will be print during the process (default
  value FALSE)}
}
\value{
  an object of class rankclust

  For example : res=rankclust(data,K=1:2,m=m)

  You can access the result by res[number of
  groups]@slotName where slotName is an element of the
  class Output.
}
\description{
  This functions runs the EM algorithm [1] or SEM-Gibbs
  algorithm for partial ranks [2] on the data in input and
  returns the estimates parameters according to ISR and
  other complementary information.
}
\examples{
data(quiz)
result=rankclust(quiz$data,K=2,m=quiz$m)
}
\author{
  Quentin Grimonprez
}
\references{
  [1] C.Biernacki and J.Jacques (2012), A generative model
  for rank data based on sorting algorithm, Computational
  Statistics and Data Analysis, 58, 162-176

  [2] J. Jacques, C. Biernacki (2012), Model-based
  clustering for multivariate partial ranking data
}
\seealso{
  See \code{\link{Output-class}} and
  \code{\link{Rankclust-class}} for available output.
}

