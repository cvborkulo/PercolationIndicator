\name{PercolationIndicator}
\alias{PercolationIndicator}
\title{
Capturing the dynamics of a network in a single parameter: the percolation indicator
}
\description{
This method captures the dynamics of a network in a single parameter: the percolation indicator. The indicator is computed by combining intensive longitudinal binary data of a system with its network structure (with temporal relationships). The dynamics of the network is modeled by two independent Poisson processes. The ratio of the parameters of these processes (infection rate and recovery rate) is used to obtain an indication of the development of the network over time. When the percolation indicator pi <= 1, activity in the network will die out. When pi > 1, the network will remain infected indefinitely. Can deal with binary data and an unweighted adjacency matrix.
}
\usage{
PercolationIndicator(data, adjacency, separate = FALSE, ...)
}

\arguments{
  \item{data}{
A matrix with multiple observations of one system. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables on a specific time point. Must be longitudinal binary data.
}
\item{adjacency}{
An unweighted adjacency matrix. The dimension of the matrix is nvars x nvars. The adjacency matrix may also be (2 x nvars) x (2 x nvars), containing temporal relationships.
}
  \item{separate}{
Logical. Can be TRUE or FALSE. The default is FALSE, resulting in one percolation indicator that applies to the whole network. When TRUE, percolation indicators per node of the network are returned.
}
  \item{\dots}{
Currently not used.
}
}

\value{
PercolationIndicator returns a "PercolationIndicator" object that contains the following items:
\item{inf.rate }{The infection rate (of the network or per node of the network).}
\item{recov.rate }{The recovery rate.}
\item{perc.ind }{The percolation indicator (of the network or per node of the network).}
}
\references{
  van Borkulo, C. D., Kamphuis, J. H., and Waldorp, L. J. (manuscript submitted for publication). Predicting behaviour of networks: the contact process as a model for dynamics.
  
  Fiocco, M. and van Zwet, W. R. (2004). Maximum likelihood estimation for the contact process. Lecture Notes-Monograph Series, 309 - 318.
  
  Grimmet, G. (2010). Probability on graphs: random processes on graphs and lattices, 1. Cambridge University Press.
  
  Harris, T. E. (1974). Contact interactions on a lattice. The Annals of Probability, 2, 969 - 988.
}
\author{
Claudia D. van Borkulo

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
\note{
See also my website: http://cvborkulo.com
}

% \examples{
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
