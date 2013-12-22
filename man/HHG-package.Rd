\name{HHG-package}
\alias{HHG-package}
\docType{package}
\title{Heller-Heller-Gorfine Tests of Independence}
\description{Heller-Heller-Gorfine (HHG) tests are a set of statistical tests of independnece between two random vectors of arbitrary dimensions given a finite sample. For testing independence between two scalar random variables, the package also contains implementations of the data-driven partitions (DDP) and all data partitions (ADP) tests which are distribution-free and thus are much faster to apply.}
\details{
\tabular{ll}{
Package: \tab HHG\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-12-22\cr
License: \tab GPL-2\cr
}
The package contains two major functions: hhg.test and xdp.test.
}
\author{
Shachar Kaufman, based on an earlier implementation of the original HHG test by Ruth Heller <ruheller@post.tau.ac.il> and Yair Heller <heller.yair@gmail.com>.
Maintainer: Shachar Kaufman <shachark@post.tau.ac.il>
}
\references{
Heller R., Heller Y., and Gorfine M. (2012). A consistent multivariate test of association based
on ranks of distances. \emph{arXiv:1201.3522v1}.

Kaufman S., Heller R., Heller Y., and Gorfine, M. (2013). Consistent distribution-free tests of association between univariate random variables. \emph{arXiv:1308.1559},
}
%\seealso{
%}
\examples{
# See examples in the documentation for hhg.test and xdp.test.
}
