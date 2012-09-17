\name{FRCC-package}
\alias{FRCC-package}
\alias{FRCC}
\docType{package}
\title{
 Fast Regularized Canonical Correlation Analysis
}
\description{
This package implements the Fast Regularized Canonical Correlation Analysis described in [Cruz-Cano et al., 2012]. 
The main idea of the algorithm is using the minimum risk estimators of the correlation matrices described in [Schafer and Strimmer, 2008] during the calculation of the Canonical correlation Structure. 
It can be considered an extesion of the work for two set of variables (blocks) mentioned in [Tenenhaus and Tenenhaus, 2011]
}
\details{
\tabular{ll}{
Package: \tab FRCC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-13\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The function frcc provides the canonical structure for two sets of variables X and Y. The rest of the functions help to visualize and interpret the values of the canonical structure.
}
\author{
Raul Cruz-Cano

Maintainer: Raul Cruz-Cano <raulcruz@umd.edu>
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.

Schafer, J; Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation
and implications for functional genomics. Statistical Applications in Genetics and Molecular
Biology 4:14, Article 32.

Tenenhaus, A.; Tenenhaus, M. (2011). Regularized Generalized Canonical Correlation Analysis.
Psychometrika 76:2, DOI: 10.1007/S11336-011-9206-8.
}
\keyword{ Fast Regularized Canonical Correlation Analysis }
\keyword{ Regularized Canonical Correlation Analysis }
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
}
