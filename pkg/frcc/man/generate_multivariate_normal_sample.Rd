\name{generate_multivariate_normal_sample}
\alias{generate_multivariate_normal_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
It generates a sample from a multinormal distribution function 
}
\description{
It generates a sample from a multinormal distribution function with the cross-covariance matrix described in [Cruz-Cano et al. 2012].
}
\usage{
generate_multivariate_normal_sample(p, q, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Number of desired variables in the dataset X.
}
  \item{q}{
Number of desired variables in the dataset Y.
}
  \item{n}{
 sample size desired.
}
}

\value{
A list of n sample units with the values for the variables of the datsets X and Y.
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
}
\author{
Raul Cruz-Cano
}

\examples{
# Example # 1 Multivariate Normal Data
p<-10
q<-10
n<-50
res<-generate_multivariate_normal_sample(p,q,n)
X<-res$X
Y<-res$Y
rownames(X)<-c(1:n)
colnames(X)<-c(1:p)
colnames(Y)<- c(1:q)
my_res<-frcc(X,Y)
dev.new()
plot_units(X,Y,my_res,1)
dev.new()
plot_variables(my_res,1,2,text_size=1.0)

## The function is currently defined as
function (p, q, n) 
{
    require("MASS")
    if ((p > 6) & (q > 6)) {
        Sigma_Z <- matrix(data = c(0), nrow <- (p + q), ncol <- (p + 
            q))
        for (i in 1:(p + q)) {
            for (j in 1:(q + p)) {
                if (i == j) {
                  Sigma_Z[i, j] <- 1
                }
                if (((i == (p + 1)) & (j == 1)) || ((j == (p + 
                  1)) & (i == 1))) {
                  Sigma_Z[i, j] <- 0.9
                }
                if (((i == (p + 2)) & (j == 2)) || ((j == (p + 
                  2)) & (i == 2))) {
                  Sigma_Z[i, j] <- 0.7
                }
                if (((i == (p + 3)) & (j == 3)) || ((j == (p + 
                  3)) & (i == 3))) {
                  Sigma_Z[i, j] <- 0.5
                }
                if (((i == (p + 4)) & (j == 4)) || ((j == (p + 
                  4)) & (i == 4))) {
                  Sigma_Z[i, j] <- 0.3
                }
                if (((i == (p + 5)) & (j == 5)) || ((j == (p + 
                  5)) & (i == 5))) {
                  Sigma_Z[i, j] <- 0.1
                }
            }
        }
        means <- matrix(data = c(0), nrow <- 1, ncol <- (p + 
            q))
        ndata <- mvrnorm(n, Sigma = Sigma_Z, mu = means)
        X <- ndata[, 1:p]
        Y <- ndata[, (p + 1):(p + q)]
        res <- list(X = matrix(data = c(0), nrow = n, ncol = p), 
            Y = matrix(data = c(0), nrow = n, ncol = q), Sigma_Z = matrix(data = c(0), 
                nrow = (p + q), ncol = (p + q)))
        res$X <- X
        res$Y <- Y
        res$Sigma_Z <- Sigma_Z
        res
    }
    else {
        print("Not enough variables")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
