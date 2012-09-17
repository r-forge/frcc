\name{off.diagonal.lambda}
\alias{off.diagonal.lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the value of the shrinkage coefficient for the off-diagonal matrices.
}
\description{
Calculates the value of the shrinkage coefficient for the off-diagonal matrices as decribed in [Cruz-Cano et al., 2012]
}
\usage{
off.diagonal.lambda(xs, p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xs}{
Matrix with the values for the datasets X and Y.
}
  \item{p}{
Number of variables in the dataset X.
}
  \item{q}{
Number of variables in the dataset Y.
}
}

\value{
Shrinkage coefficient for the off-diagonal matrices used to calculate the FRCC canonical structure.
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
}
\author{
Raul Cruz-Cano
}

\examples{
##   This is an internal function. No need for examples.

## The function is currently defined as
function (xs, p, q) 
{
    xs <- as.matrix(xs)
    w = matrix(data = c(1), ncol = dim(xs)[1], nrow = 1)
    w2 = sum(w * w)
    h1w2 = w2/(1 - w2)
    sw = sqrt(w)
    Q1.squared = (crossprod(sweep(xs, MARGIN = 1, STATS = sw, 
        FUN = "*")))^2
    Q1.squared_XY = Q1.squared[1:p, p + 1:q]
    Q2 = crossprod(sweep(xs^2, MARGIN = 1, STATS = sw, FUN = "*")) - 
        Q1.squared
    Q2_XY = Q2[1:p, p + 1:q]
    denominator = sum(Q1.squared_XY)
    numerator = sum(Q2_XY)
    if (denominator == 0) 
        lambda = 1
    else lambda = min(1, numerator/denominator * h1w2)
    return(lambda)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
