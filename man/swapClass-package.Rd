\name{swapClass-package}
\alias{swapClass-package}
\alias{swapClass}
\docType{package}
\title{
SwapClass : a simple null model adapted to abundance classes data
}
\description{
A null model randomizing semi-quantitative multi-classes (or ordinal) data by swapping sub-matrices while both the row and the column marginal sums are held constant.
}
\details{
\tabular{ll}{
Package: \tab swapClass\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-06-20\cr
License: \tab GPL-3\cr
}
}
\author{
Benjamin Borgy (benjamin.borgy@gmail.com),
Remi Perronne,
Xavier Reboud
}
\references{
Borgy B. - Dynamic and assembly of weed communities: Approach by statistical modeling (2011). Ph.D. Thesis. INRA Agroicologie & Universiti de Bourgogne. DOI: 10.13140/RG.2.1.1738.1601
}
\examples{
MAT=matrix(sample(0:3,50*20,replace=TRUE,prob=c(.7,.1,.1,.1)),ncol=20)

#Calculation of the nbPerm index
nbPerm(MAT)

#Generation of null communities
NULL_MATS=nullModel(MAT)

#First null community
NULL_MATS$sim[[1]]

#Number of times that each cell has been swapped for the first null community
NULL_MATS$perms[[1]]

#nbPerm index over the observed community and the 100 null communities
plot(c(nbPerm(MAT),unlist(lapply(NULL_MATS$sim,nbPerm))),type='l',ylab="nbPerm index")

#the number of each classes per row are equal
#between observed community and the first null community
f_table = function(x) table(factor(x,levels=0:3))
all(apply(MAT,1,f_table)==apply(NULL_MATS$sim[[1]],1,f_table))

#the number of each classes per column are equal
#between observed community and the first null community
all(apply(MAT,2,f_table)==apply(NULL_MATS$sim[[1]],2,f_table))
}
