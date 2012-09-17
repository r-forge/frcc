\name{frcc}
\alias{frcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function implements the Fast Regularized Canonical Correlation Analysis
}
\description{
This function implements the Fast Regularized Canonical Correlation algorithm described in [Cruz-Cano et al., 2012]. 
The main idea of the algorithm is using the minimum risk estimators of the correlation matrices described in [Schafer and Strimmer, 2008] during the calculation of the Canonical correlation Structure. 
It can be considered an extesion of the work for two set of variables (blocks) mentioned in [Tenenhaus and Tenenhaus, 2011]
}
\usage{
frcc(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
 numeric matrix (n by p) which contains the observations on the X variables.
}
  \item{Y}{
numeric matrix (n by q) which contains the observations on the Y variables.
}
}

\value{
A list with the follwoing components of the Canonical Structure:
 \item{cor }{Canonical correlations.}
 \item{p_values }{The corresponding p-values for the each of  the canonical correlations.}
 \item{canonical_weights_X}{The canonical weights for the variables of the dataset X.}
 \item{canonical_weights_Y}{The canonical weights for the variables of the dataset Y.}
 \item{canonical_factor_loadings_X}{The interset canonical factor loadings for the variables of the dataset X.}
 \item{canonical_factor_loadings_Y}{The interset canonical factor loadings for the variables of the dataset Y.}
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.

Schafer, J; Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation
and implications for functional genomics. Statistical Applications in Genetics and Molecular
Biology 4:14, Article 32.

Tenenhaus, A.; Tenenhaus, M. (2011). Regularized Generalized Canonical Correlation Analysis.
Psychometrika 76:2, DOI: 10.1007/S11336-011-9206-8.
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
#Example #2 Soil Specification Data
 #needed to get Ecological data
data(soilspec)
list_of_units_to_be_used<-sample(1:nrow(soilspec),14)#We will only 14 soil samples
X<- soilspec[list_of_units_to_be_used,2:9]
Y<- soilspec[list_of_units_to_be_used,10:15]
colnames(X)<-c("H. pubescens", "P. bertolonii", "T. pretense", "P. sanguisorba", "R. squarrosus", "H. pilosella", "B. media","T. drucei")
colnames(Y)<- c("d","P","K","d x P", "d x K","P x K")
my_res<-frcc(X,Y)
dev.new()
plot_units(X,Y,my_res,1) 
dev.new()
plot_variables(my_res,1,2)
#Example #3 NCI-60 micrRNA Data
data("Topoisomerase_II_Inhibitors")
data("microRNA")
my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
dev.new()
plot_variables(my_res,1,3,inner_circle_radius=.5,text_size=.65)
for( i in 1:dim(microRNA)[2])
{
 colnames(microRNA)[i]<-substr(colnames(microRNA)[i], 1, 2)
}#end for i
dev.new()
plot_units(t(microRNA),-1*t(Topoisomerase_II_Inhibitors),my_res,1,1,text_size=0.01)

## The function is currently defined as
function (X, Y) 
{
    require("corpcor")
    require("CCP")
    p <- dim(X)[2]
    q <- dim(Y)[2]
    n <- dim(X)[1]
    res.frcc <- list(cor = matrix(data = c(0), nrow = 1, ncol = q), 
        p_values = matrix(data = c(0), nrow = 1, ncol = q), canonical_weights_X = matrix(data = c(0), 
            nrow = p, ncol = q), canonical_weights_Y = matrix(data = c(0), 
            nrow = q, ncol = q))
    for (i in 1:p) {
        X[, i] <- (X[, i] - mean(X[, i]))/as.numeric(sapply(as.data.frame(X[, 
            i]), sd))
    }
    for (i in 1:q) {
        Y[, i] <- (Y[, i] - mean(Y[, i]))/as.numeric(sapply(as.data.frame(Y[, 
            i]), sd))
    }
    S_star_XX <- cor.shrink(X)
    S_star_YY <- cor.shrink(Y)
    lambda_XY <- off.diagonal.lambda(cbind(X, Y), p, q)
    S_star_XY <- cor(X, Y) * (1 - lambda_XY)
    S_star_YX <- t(S_star_XY)
    print(paste("Off-Diagonal Shrinkage Coefficient:", lambda_XY))
    latent_roots <- eigen(solve(S_star_YY) \%*\% S_star_YX \%*\% 
        solve(S_star_XX) \%*\% S_star_XY, only.values = FALSE)
    res.frcc$canonical_weights_Y <- latent_roots$vectors
    for (i in 1:q) {
        res.frcc$canonical_weights_X[, i] <- solve(S_star_XX) \%*\% 
            S_star_XY \%*\% res.frcc$canonical_weights_Y[, i]
    }
    rownames(res.frcc$canonical_weights_X) <- colnames(X)
    rownames(res.frcc$canonical_weights_Y) <- colnames(Y)
    for (i in 1:q) {
        U <- as.matrix(X) \%*\% as.numeric(res.frcc$canonical_weights_X[, 
            i])
        U <- (U - mean(U))/as.numeric(sapply(as.data.frame(U), 
            sd))
        V <- as.matrix(Y) \%*\% res.frcc$canonical_weights_Y[, 
            i]
        V <- (V - mean(V))/as.numeric(sapply(as.data.frame(V), 
            sd))
        res.frcc$cor[i] <- cor(U, V)
    }
    flag_anomaly_in_correlations <- 0
    for (i in 1:(q - 1)) {
        if (res.frcc$cor[i] < res.frcc$cor[i + 1]) {
            flag_anomaly_in_correlations <- 1
        }
    }
    if (flag_anomaly_in_correlations > 0.5) {
        print("Warning: Too little data, results might be unreliable")
        res.frcc <- rearrange.frcc(res.frcc)
    }
    p_values_struct = 0
    p_values_struct <- p.asym(res.frcc$cor, n, p, q, tstat = "Wilks")
    res.frcc$p_values <- p_values_struct$p.value
    canonical_factor_loadings_X <- matrix(data = c(0), nrow = p, 
        ncol = q)
    canonical_factor_loadings_Y <- matrix(data = c(0), nrow = q, 
        ncol = q)
    for (i in 1:q) {
        canonical_factor_loadings_X[, i] <- cor(as.matrix(X), 
            as.matrix(Y) \%*\% as.matrix(res.frcc$canonical_weights_Y[, 
                i]))
        canonical_factor_loadings_Y[, i] <- cor(as.matrix(Y), 
            as.matrix(X) \%*\% as.matrix(res.frcc$canonical_weights_X[, 
                i]))
    }
    res.frcc$canonical_factor_loadings_X <- matrix(data = c(0), 
        nrow = p, ncol = q)
    res.frcc$canonical_factor_loadings_Y <- matrix(data = c(0), 
        nrow = q, ncol = q)
    for (i in 1:q) {
        res.frcc$canonical_factor_loadings_X[, i] <- canonical_factor_loadings_X[, 
            i]
        res.frcc$canonical_factor_loadings_Y[, i] <- canonical_factor_loadings_Y[, 
            i]
    }
    rownames(res.frcc$canonical_factor_loadings_X) <- colnames(X)
    rownames(res.frcc$canonical_factor_loadings_Y) <- colnames(Y)
    print(res.frcc)
    res.frcc
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
