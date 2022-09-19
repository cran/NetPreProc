\name{Sparsify.matrix.fixed.neighbours-methods}
\docType{methods}
%\alias{Sparsify.matrix.fixed.neighbours-methods}
\alias{Sparsify.matrix.fixed.neighbours,graph-method}
\alias{Sparsify.matrix.fixed.neighbours,matrix-method}
\alias{Sparsify.matrix.fixed.neighbours}
\title{Sparsifying the graph by a fixed number of edges per node}
\description{
 Methods to sparsify a network matrix by fixing the number of edges for each node.
It selects the first k neighbours for each node (by row) according to the weight of the edge
By this function you select exactly k edges for each node (if there are at least k edges in the adjacency matrix).
The resulting matrix is not symmetric.
}

\usage{
Sparsify.matrix.fixed.neighbours(W, k=10)
}
\arguments{
  \item{W}{
 an object representing the graph to be normalized
}
  \item{k}{
the number of  edges for each node (def.=10)
}
}
\value{
a sparsified matrix (Warning: the matrix is not symmetric)
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
# sparsification with 10 neighbours per node 
W10 <- Sparsify.matrix.fixed.neighbours(W);
# sparsification with 20 neighbours per node 
\donttest{W20 <- Sparsify.matrix.fixed.neighbours(W, k=20);}
}


\keyword{methods}
\keyword{graph pre-processing}
