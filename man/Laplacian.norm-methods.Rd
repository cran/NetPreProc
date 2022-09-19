\name{Laplacian.norm-methods}
\docType{methods}
%\alias{Laplacian.norm-methods}
\alias{Laplacian.norm,graph-method}
\alias{Laplacian.norm,matrix-method}
\alias{Laplacian.norm}
\title{ Laplacian graph normalization }
\description{
 Methods to normalize weights of square symmetric adjacency matrices.
 A network matrix is normalized by dividing each entry \eqn{W_{ij}} by the square root of the product of the sum of elements of row \eqn{i} and the sum of the elemnts in column \eqn{j}.
In other words if \eqn{D} is a diagonal matrix such that \eqn{D_{ii} = \sum_j W_{ij}}, then the normalize matrix is: 
\deqn{
W_{norm} = D^{-1/2}  W  D^{-1/2}
}
}

\usage{
Laplacian.norm(W)
%%\S4method{Binary.matrix.by.thresh}{"matrix"}(W, thresh=0.5)
}
\arguments{
  \item{W}{
 an object representing the graph to be normalized
}
}
\value{
The normalized adjacency matrix of the network
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
# normalization of drug-drug similarity networks
data(DD.chem.data);
W <- Laplacian.norm(DD.chem.data);
# the same using an object of class graphAM
\donttest{g <- new("graphAM", adjMat=DD.chem.data, values=list(weight=DD.chem.data));
Wg <- Laplacian.norm(g);}
}


\keyword{methods}
\keyword{graph normalization}
