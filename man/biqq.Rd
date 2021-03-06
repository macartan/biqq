% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biqq.R
\name{biqq}
\alias{biqq}
\title{biqq function, taking data as arguments
Basic model; used to define fast biqq function}
\usage{
biqq(fit = fit_biqq, XY = c(0, 0, 0, 0), XYK = c(0, 0, 0, 0, 0, 0, 0,
  0), XYK_data = NULL, alpha_prior = c(1, 1, 1, 1), pi_alpha = c(1,
  1, 1, 1), pi_beta = c(1, 1, 1, 1), q0_alpha = c(1, 1, 1, 1),
  q0_beta = c(1, 1, 1, 1), q1_alpha = c(1, 1, 1, 1), q1_beta = c(1,
  1, 1, 1), doubly_decisive = FALSE, iter = 1000, chains = 2, ...)
}
\arguments{
\item{fit}{a fitted model}

\item{XY}{XY data -- 4-vector  00 01 10 11}

\item{XYK}{XYK data -- 12-vector 000 001 010 011...}

\item{XYK_data}{XYK data provided as a dataframe}

\item{doubly_decisive}{if TRUE pobative value defaults to clue traces for a and b types only}
}
\value{
A stan object
}
\description{
biqq function, taking data as arguments
Basic model; used to define fast biqq function
}
\examples{
print(biqq(alpha_prior = c(5,1,1,1)))

}
\keyword{biqq}
