\name{Prob.norm-methods}
\docType{methods}
%\alias{Prob.norm-methods}
\alias{Prob.norm,graph-method}
\alias{Prob.norm,matrix-method}
\alias{Prob.norm}
\title{ Probabilistic normalization of a graph }
\description{
Method to compute the transition probability matrix of network.
A network matrix is normalized by dividing each entry \eqn{W_{ij}} by the the sum of elements of row \eqn{i} 
In other words if \eqn{D} is a diagonal matrix such that \eqn{D_{ii} = \sum_j W_{ij}} then the normalize matrix is: 
\deqn{
W_{norm} = D^{-1}  W 
}
}

\usage{
Prob.norm(W)
%%\S4method{Binary.matrix.by.thresh}{"matrix"}(W, thresh=0.5)
}
\arguments{
  \item{W}{
 an object representing the graph to be normalized
}
}
\value{
The normalized transition probability matrix of network
}

\section{Methods}{
\describe{

\item{\code{signature(W = "graph")}}{
an object of the virtual class graph (hence including objects of class \code{\link[graph:graphAM-class]{graphAM}}  and \code{\link[graph:graphNEL-class]{graphNEL}} from the package \pkg{graph})
}

\item{\code{signature(W = "matrix")}}{
a matrix representing the adjacency matrix of the graph
}
}}


\examples{
library(bionetdata);
# making transition prob matrix for a drug-drug similarity network
data(DD.chem.data);
W <- Prob.norm(DD.chem.data);
# the same using an object of class graphAM and of class graphNEL
\donttest{
g <- new("graphAM", adjMat=DD.chem.data, values=list(weight=DD.chem.data));
Wg <- Prob.norm(g);
g2 <- as(g, "graphNEL");
Wg2 <- Prob.norm(g2);}
}


\keyword{methods}
\keyword{graph pre-processing}
