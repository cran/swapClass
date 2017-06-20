\name{nbPerm}
\alias{nbPerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Number of swaps that can be performed
}
\description{
The number of permutations is calculated by summing the number of "swappable" sub-matrices per pair of abundance classes. The number of abundance classes needs to be relatively low, the swapClass algorithm being efficient on semi-quantitative multi-classes (or ordinal) variables (e.g. Braun-Blanquet- type abundance/dominance scores) but not on continuous quantitative variables.
}
\usage{
nbPerm(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
an object of class matrix with abundance classes.
}
}
\details{
nbPerm is the sum of "swappable" sub-matrices divided by the number of cells in the matrix (number of rows x number of columns). 
nbPerm calls a dynamically defined functions (swapClass_NBP) with in-lined C code compiled by \link[inline]{setCMethod}  on package load.
}
\value{
a single numeric value.
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
\code{\link{nullModel}}
}
\examples{
MAT=matrix(sample(0:3,50*20,replace=TRUE,prob=c(.7,.1,.1,.1)),ncol=20)
nbPerm(MAT)
}

