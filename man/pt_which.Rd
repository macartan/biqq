% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_tracing_functions.R
\name{pt_which}
\alias{pt_which}
\title{A function to figure out which worlds satisfy some query}
\usage{
pt_which(model, U = NULL, ops, query, sims = 500)
}
\arguments{
\item{model}{A model made using \code{\link{pt_model}}.}

\item{U}{A matrix of contexts (worlds). If missing this is generated based on sims.}

\item{ops}{A set of operations.}

\item{query}{A query defined on outcomes on observables, V.}

\item{sims}{Optional number of simulations for draws; defaults to 500 if U not defined.}
}
\description{
The output is a list that contains a matrix of contexts, U, the outcomes associated with them, V, and an indicator for whether the query is satisfied or not in each context, A.
}
\examples{
pt_which(model = pt_model(),
         ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA)),
         query = function(x) (x[[1]][5] == 1) & (x[[2]][5] == 0),
         sims = 3)

}
