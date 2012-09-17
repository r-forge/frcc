\name{plot_units}
\alias{plot_units}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the experimental units in the Canonical Variates Space
}
\description{
 This function plots the experimental units used in the FRCCA as points in a two-dimensional 
 plane in which the axis are the canonical variates selected by the user
}
\usage{
plot_units(X, Y, res.mrcc, i, text_size = 0.8, point_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
numeric matrix (n by p) which contains the observations on the X variables.
}
  \item{Y}{
numeric matrix (n by p) which contains the observations on the Y variables.
}
  \item{res.mrcc}{
List containing a canonical structure provided by the function frcc for the dataset X and Y.
}
  \item{i}{
Canonical Variate which will be used for the axes (X for horizontal and Y for vertical).
}
  \item{text_size}{
	Character expansion factor for the labels of the experimental units.
}
  \item{point_size}{
	Character expansion factor for the point representing the experimental units.
}
}

\value{

This function just creates the units plot. It does not return a value.
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
#Example #3 NCI-60 micrRNA Data
data("Topoisomerase_II_Inhibitors")
data("microRNA")
my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
for( i in 1:dim(microRNA)[2])
{
 colnames(microRNA)[i]<-substr(colnames(microRNA)[i], 1, 2)
}#end for i
dev.new()
plot_units(t(microRNA),-1*t(Topoisomerase_II_Inhibitors),my_res,1,1,text_size=0.01)

## The function is currently defined as
function (X, Y, res.mrcc, i, text_size = 0.8, point_size = 2) 
{
    require("calibrate")
    U <- as.matrix(X) \%*\% as.numeric(res.mrcc$canonical_weights_X[, 
        i])
    U <- (U - mean(U))/as.numeric(sapply(as.data.frame(U), sd))
    V <- as.matrix(Y) \%*\% res.mrcc$canonical_weights_Y[, i]
    V <- (V - mean(V))/as.numeric(sapply(as.data.frame(V), sd))
    s1 <- paste("U", i, sep = "")
    s2 <- paste("r=", round(res.mrcc$cor[i], digits = 3), sep = "")
    s3 <- paste("p-value=", round(res.mrcc$p_value[i], digits = 3))
    my_xlab = paste(s1, s2, s3, sep = "       ")
    my_ylab = paste("V", i)
    plot(U, V, ylim = c(min(V) - 0.5, max(V) + 0.5), xlim = c(min(U) - 
        0.5, max(U) + 0.5), col = "white", xlab = my_xlab, ylab = my_ylab)
    points(U, V, cex = point_size, pch = 19, col = "black")
    textxy(U, V, rownames(X), cx = text_size, dcol = "black")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
