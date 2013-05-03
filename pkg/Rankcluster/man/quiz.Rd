\docType{data}
\name{quiz}
\alias{quiz}
\title{Multidimensionnal rank data : quiz}
\format{A list containing:

-data: A matrix of size 70*16. Each row correspond to the answer of a student, a answer corresponds to 4 ranks (ordering representation) of size 4.

-frequency: A matrix of size 63*17. Each rows corresponds to one of the 63 differents observed rank (ordering representation). A row contains 4 ranks of size 4 and the frequency.

-m: a vector with the size of the rank for each dimension.}
\source{
  Creator: Julien Jacques
}
\description{
  70 students (40 of third year and 30 of fourth year of statistics engineering) from Polytech'Lille have to answer 4 quizzes:


1) Literature Quiz
This quiz consists of ranking four french writers according to chronological order:
 1=Victor Hugo, 2=Moliere, 3=Albert Camus, 4=Jean-Jacques Rousseau.

2) Football Quiz
This quiz consists of ranking four national football teams according to increasing number of wins in thefootball World Cup: 1=France, 2=Germany, 3=Brazil, 4=Italy.

3) Mathematics Quiz
This quiz consists of ranking four numbers according to increasing order: 
1=pi/3, 2=log(1), 3=exp(2), 4=(1+sqrt(5))/2.

4) Cinema Quiz
This quiz consists of ranking four Tarentino's movies according to chronological order: 
1=Inglourious Basterds, 2=Pulp Fiction, 3=Reservoir Dogs, 4=Jackie Brown.
	
}
\examples{
data(quiz)
}
\keyword{datasets}

