\docType{data}
\name{APA}
\alias{APA}
\title{rank data : APA}
\format{A list containing :

-data: matrix of size 5738x5 containing the 5738 full observed ranks in orderinfg representation. The ordering representation x = (x_1,...,x_m ) means that Object O_x_i is in the ith position (i =1,...,m); For example, if the ordering rank is (4,3,1,2,5), it means that judge ranks first tthe cadidate 4, second  the candidate 3, ....


-frequency: matrix of size 120x6. Each row corresponds to one of the different observed rank. The first fifth columns contains the observed ranks (ordering representation) and the sisth column contains the frequency of observation.

-m: vector with the size of the rank for each dimension.}
\source{
  "Group representations in probability and statistics",P. Diaconis,1988
}
\description{
This data contains 5738 votes of member of American Psychological Association (APA) presidential election of 1980.
Membres had to rank five candidates in order of preference. There is a total of 15449 votes with 5738 full ranking).

For information, candidates 1 and 3 were research psychologists, candidates 4 and 5 were clinical psychologists and candidate 2 was a community psychologist.
	
}
\examples{
data(APA)
}
\keyword{datasets}

