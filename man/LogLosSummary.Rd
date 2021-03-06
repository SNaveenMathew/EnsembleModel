% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caretEnsemble_package.R
\name{LogLosSummary}
\alias{LogLosSummary}
\title{Binary log loss}
\usage{
LogLosSummary(data, lev = NULL, model = NULL)
}
\arguments{
\item{data}{Vector of observed classes}

\item{lev}{(Optional) Levels of the outcome variable}

\item{model}{Model predictions for p(y=1)}
}
\value{
List of accuracy, kappa and binary logloss
}
\description{
This function calculates binary logloss given the actual classes and
 predictions
}
\examples{
LogLossSummary(data = c(0,1), model = c(0.1, 0.9))

}
