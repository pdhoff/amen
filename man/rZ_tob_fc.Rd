% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rZ_tob_fc.R
\name{rZ_tob_fc}
\alias{rZ_tob_fc}
\title{Simulate Z based on a tobit model}
\usage{
rZ_tob_fc(Z, EZ,rho,s2,Y)
}
\arguments{
\item{Z}{a square matrix, the current value of Z}

\item{EZ}{expected value of Z}

\item{rho}{dyadic correlation}

\item{s2}{dyadic variance}

\item{Y}{square relational matrix with nonnegative entries}
}
\value{
a square matrix, the new value of Z
}
\description{
Simulates a random latent matrix Z given its expectation, dyadic correlation
and a nonnegative relational matrix Y
}
\author{
Peter Hoff
}
