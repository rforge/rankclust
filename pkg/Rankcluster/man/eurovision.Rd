\docType{data}
\name{eurovision}
\alias{eurovision}
\title{Multidimensionnal partial rank data : eurovision}
\format{A list containig:

-data: A matrix of size 34*48. A row contains a multidimensionnal rank in ordering representation. The first eigthth rows corresponds to the 2007 contest, the next eigth to the 2008 contest, etc...

-frequency: A matrix of size 34*49 containing the differents multidimensionnal ranks. A row contains 6 ranks of size 8 and the last column contains the frequency

-m: a vector with the size of ranks for each dimension.}
\source{
  http://www.eurovision.tv
}
\description{
The data contains the ranking by 33 countries of 8 countries during Eurovision song contest from 2007 to 2012 by decreasing number of points given.


The 8 countries are the common finalists from 2007 to 2012 constest, there are:

1: France, 2:Germany, 3:Greece, 4:Romania, 5:Russia, 6:Spain, 7:Ukrain, 8:United Kingdom.
They are ranking according to the number of points given by a voter. 

The 33 voters are the common voters countries to these six contests. 

If a voter has ranked only 4 countries then the rank of size m contains the number of these 4 countries and four zeros.
	
}
\examples{
data(eurovision)
}
\keyword{datasets}

