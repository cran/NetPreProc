\name{Magnify.binary.features.norm-methods}
\docType{methods}
%\alias{Magnify.binary.features.norm-methods}
\alias{Magnify.binary.features.norm,matrix-method}
\alias{Magnify.binary.features.norm}
\title{ Normalization of binary matrices }
\description{
 Normalization of binary matrices according to the procedure described in Mostafavi et al. 2008.
Having a binary matrix M, for each feature, if b is the proportion of 1, then ones are replaced with -log(b) and zeros with log(1-b).
}

\usage{
Magnify.binary.features.norm(M)
%%\S4method{Binary.matrix.by.thresh}{"matrix"}(W, thresh=0.5)
}
\arguments{
  \item{M}{
 an object representing the matrix to be normalized
}
}
\value{
The normalized  matrix 
}

\section{Methods}{
\describe{

\item{\code{signature(M = "matrix")}}{
Input binary matrix. Rows are examples, columns features
}
}}


\examples{
D <- matrix(ifelse(runif(40000)>0.9,1,0),nrow=100);
M <- Magnify.binary.features.norm(D);
}

\references{
Mostafavi, S., Ray, D., Warde-Farley, D., Grouios, C., & Morris, Q.
GeneMANIA: a real-time multiple association network integration algorithm
for predicting gene function. Genome Biology, 9, 2008.
}


\keyword{methods}
