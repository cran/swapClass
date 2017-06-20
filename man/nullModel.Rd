\name{nullModel}
\alias{nullModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SwapClass null model
}
\description{
 a simple null model adapted to semi-quantitative multi-classes (or partially ordered) data.
}
\usage{
nullModel(mat, nsim = 100, recursive = TRUE, burnin = NULL, thin = NULL, strata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
an object of class matrix (or data.frame) with numerical abundance classes (site by row and species by column).
}
  \item{nsim}{
Number of simulated null communities (positive integer), default = 100.
}
  \item{recursive}{
A logical value indicating if generation of a new null community is performed from the last one created, defaut = T.
}
  \item{burnin}{
Number of null communities discarded before proper analysis (positive integer), if NULL (defaut) then burnin = nrow(mat)*ncol(mat)*10.
}
  \item{thin}{
Number of discarded null communities between two evaluations (positive integer) if recursive=T, if NULL (defaut) then thin = max(nrow(mat)*ncol(mat)*10,1000).
}
  \item{strata}{
a numeric vector of length egal to nrow(mat) supplying strata. Swap are performed within strata, defaut = NULL.
}
}
\details{
 The SwapClass model is derived from the "swap philosophy" apply on presence-absence data. 
The limited number of abundance classes and their repetitiveness allow the extension of the "swap" method classically applied on presence-absence data with two classes (0 and 1). In the same way that "swap" methods permute sub-matrices of presence/absence community matrix, semi-quantitative data can be randomized by swapping sub-matrices while row and column marginals are not modified.
nullModel calls a dynamically defined functions (swapClass_swapC) with in-lined C code compiled by \link[inline]{setCMethod} on package load.
}
\value{
\item{sim}{an object of class list containing the null communities}
\item{perms}{an object of class list containing matrices with the number of times each cell has been swapped.}
}
\author{
Benjamin Borgy (benjamin.borgy@gmail.com),
Remi Perronne,
Xavier Reboud
}
\references{
Borgy B. - Dynamic and assembly of weed communities: Approach by statistical modeling (2011). Ph.D. Thesis. INRA Agroicologie & Universiti de Bourgogne. DOI: 10.13140/RG.2.1.1738.1601
}
\seealso{
\code{\link{nbPerm}}
}
\examples{
MAT=matrix(sample(0:3,50*20,replace=TRUE,prob=c(.7,.1,.1,.1)),ncol=20)

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

