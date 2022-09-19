\name{check.network-methods}
\docType{methods}
%\alias{check.network-methods}
\alias{check.network,graph-method}
\alias{check.network,matrix-method}
\alias{check.network}
\title{ Graph checking }
\description{
 Method to check the characteristics of a graph.
Check if its adjacency matrix is symmetric, if it has NA, NaN o Inf values, and some minimals 
statistics about nodes and edges.
}

\usage{
check.network(W, name="Network matrix")
}
\arguments{
  \item{W}{
 an object representing the graph to be checked
}
  \item{name}{
a character vector that will be printed as heading
}
}
\value{
It return a list of strings about the characteristics of the graph
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
data(DD.chem.data);
check.network(DD.chem.data);
W <- Prob.norm(DD.chem.data);
check.network(W, "prob. transition matrix");
\donttest{WL <- Laplacian.norm(DD.chem.data);
check.network(WL, "Laplacian norm. matrix");}

library(graph)
g1 = randomEGraph(LETTERS[1:15], edges = 40);
check.network(g1, "random graph");
}


\keyword{methods}
\keyword{utilities}
