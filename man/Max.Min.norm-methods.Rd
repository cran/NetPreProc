\name{Max.Min.norm-methods}
\docType{methods}
%\alias{Max.Min.norm-methods}
\alias{Max.Min.norm,graph-method}
\alias{Max.Min.norm,matrix-method}
\alias{Max.Min.norm}
\title{ Max-min graph normalization }
\description{
Graph normalization with respect to the minimum and maximum value of its weights.
Each entry of the normalized matrix is in the range \eqn{[0..1]}:

\deqn{W_{norm} = \frac{(W - \min(W))}{(\max(W)-\min(W))}}{W.norm = (W - min(W))/(max(W)-min(W))} 
}

\usage{
Max.Min.norm(W)
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
# max-min normalization for a drug-drug similarity network
data(DD.chem.data);
W <- Max.Min.norm(DD.chem.data);
# the same using an object of class graphAM 
\donttest{
g <- new("graphAM", adjMat=DD.chem.data, values=list(weight=DD.chem.data));
Wg <- Max.Min.norm(g);}
}


\keyword{methods}
\keyword{graph normalization}
