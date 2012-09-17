\name{plot_variables}
\alias{plot_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variables in the Canonical Factor Loadings Space
}
\description{
 This function plots the variables used in the FRCCA as points in a two-dimensional 
 plane in which the axis are the canonical factor loadings selected by the user.
}
\usage{
plot_variables(res.mrcc, i, j, inner_circle_radius = 0.5, text_size = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.mrcc}{
  List containing a canonical structure provided by the function frcc.
}
  \item{i}{
 Canonical Factor Loadings which will be used as the horizontal axis.
}
  \item{j}{
  Canonical Factor Loadings  which will be used as the vertical axis.
}
  \item{inner_circle_radius}{
  Radius of the circle which is used to determine which variables are significant. Only the significant variables will be labled. 
}
  \item{text_size}{
 	Character expansion factor for the labels of the variables.
}
}

\value{
 This function just creates the variables plot. It does not return a value.
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
plot_variables(my_res,1,2)
#Example #3 NCI-60 micrRNA Data
data("Topoisomerase_II_Inhibitors")
data("microRNA")
my_res <- frcc(t(microRNA),t(-1*Topoisomerase_II_Inhibitors))
dev.new()
plot_variables(my_res,1,3,inner_circle_radius=.5,text_size=.65)

## The function is currently defined as
function (res.mrcc, i, j, inner_circle_radius = 0.5, text_size = 0.8) 
{
    require("calibrate")
    p <- nrow(res.mrcc$canonical_factor_loadings_X)
    q <- nrow(res.mrcc$canonical_factor_loadings_Y)
    first_dimension <- matrix(data = c(0), nrow = p + q, ncol = 1)
    first_dimension[1:p, 1] <- res.mrcc$canonical_factor_loadings_X[, 
        i]
    first_dimension[(p + 1):(p + q), 1] <- res.mrcc$canonical_factor_loadings_Y[, 
        i]
    second_dimension <- matrix(data = c(0), nrow = p + q, ncol = 1)
    second_dimension[1:p, 1] <- res.mrcc$canonical_factor_loadings_X[, 
        j]
    second_dimension[(p + 1):(p + q), 1] <- res.mrcc$canonical_factor_loadings_Y[, 
        j]
    my_xlab = paste("CC", i)
    my_ylab = paste("CC", j)
    plot(first_dimension, second_dimension, ylim = c(-1.25, 1.25), 
        xlim = c(-1.25, 1.25), pch = " ", xlab = my_xlab, ylab = my_ylab)
    points(first_dimension[1:p, ], second_dimension[1:p, ], cex = 2, 
        pch = 17, col = "grey70")
    points(first_dimension[(p + 1):(p + q), ], second_dimension[(p + 
        1):(p + q), ], cex = 2, pch = 19, col = 1)
    points_names <- c(rownames(res.mrcc$canonical_factor_loadings_X), 
        rownames(res.mrcc$canonical_factor_loadings_Y))
    for (i in 1:(p + q)) {
        if (sqrt((first_dimension[i, 1]^2) + (second_dimension[i, 
            1]^2)) < inner_circle_radius) {
            points_names[i] <- " "
        }
    }
    textxy(first_dimension, second_dimension, points_names, cx = text_size, 
        dcol = "blue")
    custom.draw.circle(0, 0, inner_circle_radius, col = "black")
    custom.draw.circle(0, 0, 1, col = "black")
    segments(-1.45, 0, 1.45, 0)
    segments(0, -1.45, 0, 1.45)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
