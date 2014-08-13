\name{PercolationIndicator-package}
\alias{PercolationIndicator-package}
\docType{package}
\title{
  Capturing the dynamics of a network in a single parameter: the percolation indicator
}
\description{
  This method captures the dynamics of a network in a single parameter: the percolation indicator. The indicator is computed by combining intensive longitudinal binary data of a system with its network structure (with temporal relationships). The dynamics of the network is modeled by two independent Poisson processes. The ratio of the parameters of these processes (infection rate and recovery rate) is used to obtain an indication of the development of the network over time. When the percolation indicator pi <= 1, activity in the network will die out. When pi > 1, the network will remain infected indefinitely. Can deal with binary data and an unweighted adjacency matrix.
}
\details{
  \tabular{ll}{
    Package: \tab PercolationIndicator\cr
    Type: \tab Package\cr
    Version: \tab 0.1.0\cr
    Date: \tab 2014-12-08\cr
    License: \tab What license is it under?\cr
  }
}
\author{
  Claudia D. van Borkulo
  
  Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
\references{
  van Borkulo, C. D., Kamphuis, J. H., and Waldorp, L. J. (manuscript submitted for publication). Predicting behaviour of networks: the contact process as a model for dynamics.
  
  Fiocco, M. and van Zwet, W. R. (2004). Maximum likelihood estimation for the contact process. Lecture Notes-Monograph Series, 309 - 318.
  
  Grimmet, G. (2010). Probability on graphs: random processes on graphs and lattices, 1. Cambridge University Press.
  
  Harris, T. E. (1974). Contact interactions on a lattice. The Annals of Probability, 2, 969 - 988.
}
% ~~ Optionally other standard keywords, one per line, from file ~~
  % ~~ KEYWORDS in the R documentation directory ~~