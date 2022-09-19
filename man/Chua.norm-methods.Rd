\name{Chua.norm-methods}
\docType{methods}
%\alias{Chua.norm-methods}
\alias{Chua.norm,graph-method}
\alias{Chua.norm,matrix-method}
\alias{Chua.norm}
\title{ Chua normalization }
\description{
 Normalization of graphs according to Chua et al., 2007. 
 The normalized weigths  between nodes are computed by taking into account their neighborhoods.
 This normalization is meaningful in particular with interaction data.
 More precisely, the normalized weigth \eqn{W_{ij}} between nodes \eqn{i} and \eqn{j} is computed by taking into account their neighborhods \eqn{N_i} and \eqn{N_j} :

\deqn{W_{ij} = \frac{2|N_i \cap N_j|}{|N_i \setminus N_j| + 2|N_i \cap N_j| + 1}\times \frac{2|N_i \cap N_j|}{|N_j \setminus N_i| + 2|N_i \cap N_j| + 1}}
where \eqn{N_k} is the set of the neighbors of gene \eqn{k} (\eqn{k} is included).   
}

\usage{
Chua.norm(W)
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
\donttest{library(bionetdata);
data(Yeast.Biogrid.data);
W <- Chua.norm(Yeast.Biogrid.data);}
}

\references{
Chua, H., Sung, W., & Wong, L. An efficient strategy for extensive
integration of diverse biological data for protein function prediction. Bioinformatics, 23 , 3364--3373, 2007.
}

\keyword{methods}
\keyword{graph normalization}
