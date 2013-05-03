\docType{data}
\name{words}
\alias{words}
\title{rank data : words}
\format{A list containing :

-data: A matrix of size 98*5 containig the 98 answers. each row is a rank (ordering representation)

-frequency: matrix of size 15*6. The matrix contains the 15 differents observed ranks. A row contains 1 rank of size 5 and the frequency.

-m: vector with the size of the rank for each dimension.}
\source{
M.A. Fligner and J.S. Verducci. "Distance based ranking models". J. Roy. Statist.
Soc. Ser. B, 48(3):359–369, 1986.


}
\description{
The data was collected under the auspices of the Graduate Record
Examination Board. A sample of 98 college students were asked to rank five words according to strength of association (least to most associated) with the target word "Idea":
1 = Thought, 2 = Play, 3 = Theory, 4 = Dream and 5 = Attention
	
}
\examples{
data(words)
}
\keyword{datasets}

