\name{Do.sim.matrix.Pearson}
\alias{Do.sim.matrix.Pearson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construction of the Pearson correlation matrix }
\description{
Function to obtain the Pearson correlation matrix between rows of a given matrix.
}
\usage{
Do.sim.matrix.Pearson(M, cut = TRUE, remove.negatives = TRUE, min.thresh = 0)
}
\arguments{
  \item{M}{
input matrix
}
  \item{cut}{
if TRUE (def.) at least one edge is maintained for each node, all the other edges are set to 0. If false
      no edgeis set to 0.
}
  \item{remove.negatives}{
if TRUE (def) negative values are replaced with 0 in the correlation matrix
}
  \item{min.thresh}{
minimum allowed threshold (def. 0).
If a threshold lower than min.thresh is selected, thanit is substituted by min.thresh. 
Warning: setting min.thresh to large values may lead to highly disconneted network
} 
}
\details{
You can also "sparsify" the matrix, by putting to 0 all the weights, by setting a threshold
such that at least one edge is maintained for each node.
The diagonal values are set to 0.
}
\value{
a square symmetric matrix of the Pearson correlation coefficients 
computed between the rows of M
}

\examples{
# a gaussian random matrix
D <- matrix(rnorm(20000),nrow=200);
W <- Do.sim.matrix.Pearson (D);
# the same without default parameters
W2 <- Do.sim.matrix.Pearson (D, cut=FALSE, remove.negatives=FALSE, min.thresh=-20);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graph pre-processing}

