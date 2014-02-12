\docType{methods}
\name{[}
\alias{[}
\alias{[,Rankclust-method}
\title{Getter method for rankclust output}
\usage{
\S4method{[}{Rankclust}(x, i, j, drop)
}
\arguments{
  \item{x}{object from which to extract element(s) or in
  which to replace element(s).}

  \item{i}{the number of cluster of the element we want to
  extract.}

  \item{drop}{For matrices and arrays.  If \code{TRUE} the
  result is coerced to the lowest possible dimension (see
  the examples).  This only works for extracting elements,
  not for the replacement.  See \code{\link{drop}} for
  further details.  }
}
\description{
This is overloading of square braces to extract values of
various slots of the output from the function
\code{\link{rankclust}}.
}

