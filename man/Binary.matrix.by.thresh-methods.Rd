\name{Binary.matrix.by.thresh-methods}
\docType{methods}
%\alias{Binary.matrix.by.thresh-methods}
\alias{Binary.matrix.by.thresh,graph-method}
\alias{Binary.matrix.by.thresh,matrix-method}
\alias{Binary.matrix.by.thresh}
\title{ Transforming a real-valued network matrix into a binary matrix }
\description{
Methods to transform a a real-valued network matrix into a binary matrix.
The binary matrix is obtained by thresholding: values above the given threshold are set to 1, otherwise to 0
}

\usage{
Binary.matrix.by.thresh(W, thresh=0.5)
%%\S4method{Binary.matrix.by.thresh}{"matrix"}(W, thresh=0.5)
}
\arguments{
  \item{W}{
 an object representing the graph to be normalized
}
  \item{thresh}{
the threshold (def.=0.5)
}
}
\value{
The normalized binary adjacency matrix of the network
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
W <- Binary.matrix.by.thresh(DD.chem.data);

# Using both methods with both signatures "matrix" and "graph"
# reducing dimension of the graph
library(graph);
DD.chem.data.red <- DD.chem.data[1:100,1:100];
W.red <- Binary.matrix.by.thresh(DD.chem.data.red);
g <- new("graphAM", adjMat=DD.chem.data.red, values=list(weight=DD.chem.data.red));
Wg <- Binary.matrix.by.thresh(g);
any(W.red!=Wg);
}


\keyword{methods}
\keyword{graph pre-processing}
