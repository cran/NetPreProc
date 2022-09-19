\name{Sparsify.matrix-methods}
\docType{methods}
%\alias{Sparsify.matrix-methods}
\alias{Sparsify.matrix,graph-method}
\alias{Sparsify.matrix,matrix-method}
\alias{Sparsify.matrix}
\title{Sparsifying the graph}
\description{
 Methods to sparsify a network matrix.
 By this method a general threshold is et such that you  a minimum of k edges is guranteed for each node
}

\usage{
Sparsify.matrix(W, k=1)
%%\S4method{Binary.matrix.by.thresh}{"matrix"}(W, thresh=0.5)
}
\arguments{
  \item{W}{
 an object representing the graph to be sparsified
}
  \item{k}{
the number of guaranteed edges for each node (def.=1)
}
}
\value{
The sparsified adjacency matrix of the network
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
data(FIN.data);
W <- Laplacian.norm(as.matrix(FIN.data));
# sparsification by maintaining at least one neighbour per node 
W1 <- Sparsify.matrix(W);
# sparsification by maintaining at least 20 neighbours per node (if any)
\donttest{W20 <- Sparsify.matrix(W, k=20);}

}


\keyword{methods}
\keyword{graph pre-processing}

