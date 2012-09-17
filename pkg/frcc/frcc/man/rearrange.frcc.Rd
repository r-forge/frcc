\name{rearrange.frcc}
\alias{rearrange.frcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Rearranges the canonical structure according to the canonical correlations
}
\description{
 By using the minimum risk estimators of the correlation matrices instead of the 
 sample correlation matrices the FRCC algoeithm might disrupt the order of the
 canonical correlations and hence of the canonical structure. This is unacceptable for the algorithm used  to calculate the p-values 
 which requires the canonical correltions to be ordered in a descending order. This function rearranges the canonical structure according 
 to the canonical correlations from largest to smallest. 
}
\usage{
rearrange.frcc(res.frcc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.frcc}{
  
  List containing a canonical structure produced by the function frcc.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{res.frcc }{List containing the sorted canonical structure.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review.
}
\author{
Raul Cruz-Cano
}

\examples{
##   This is an internal function. No need for examples.

## The function is currently defined as
function (res.frcc) 
{
    p <- dim(res.frcc$canonical_weights_X)[1]
    q <- dim(res.frcc$canonical_weights_Y)[1]
    aux.res.frcc <- list(cor = 0, canonical_weights_X = matrix(data = c(0), 
        nrow = p, ncol = 1), canonical_weights_Y = matrix(data = c(0), 
        nrow = q, ncol = 1))
    for (i in 1:(q - 1)) {
        for (j in (i + 1):q) {
            if (res.frcc$cor[j] > res.frcc$cor[i]) {
                aux.res.frcc$cor <- res.frcc$cor[i]
                aux.res.frcc$canonical_weights_X <- res.frcc$canonical_weights_X[, 
                  i]
                aux.res.frcc$canonical_weights_Y <- res.frcc$canonical_weights_Y[, 
                  i]
                res.frcc$cor[i] <- res.frcc$cor[j]
                res.frcc$canonical_weights_X[, i] <- res.frcc$canonical_weights_X[, 
                  j]
                res.frcc$canonical_weights_Y[, i] <- res.frcc$canonical_weights_Y[, 
                  j]
                res.frcc$cor[j] <- aux.res.frcc$cor
                res.frcc$canonical_weights_X[, j] <- aux.res.frcc$canonical_weights_X
                res.frcc$canonical_weights_Y[, j] <- aux.res.frcc$canonical_weights_Y
            }
        }
    }
    res.frcc
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
