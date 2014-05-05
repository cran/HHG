\name{HHG-package}
\alias{HHG-package}
\alias{HHG}
\docType{package}
\title{Heller-Heller-Gorfine Tests of Independence}
\description{Heller-Heller-Gorfine (HHG) tests are a set of powerful statistical tests of independnece between two random vectors of arbitrary dimensions, given a finite sample. For testing independence between two scalar random variables, the package also contains implementations of the data-derived partitions (DDP) and all-data-partitions (ADP) tests, which are distribution-free and thus much faster to apply.}
\details{
\tabular{ll}{
Package: \tab HHG\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2014-05-05\cr
License: \tab GPL-2\cr
}
The package contains two major functions: \code{hhg.test}, which implements the HHG multivariate independence test described in Heller et al. (2013), and \code{xdp.test}, which implements the DDP and ADP univariate independence tests from Heller et al. (2014). In addition, in order to demonstrate the tests, the function \code{hhg.example.datagen} generates simulated data, and the dataset \code{hughes} contains some gene expression data.
}
\author{
Shachar Kaufman, based in part on an earlier implementation of the original HHG test by Ruth Heller <ruheller@post.tau.ac.il> and Yair Heller <heller.yair@gmail.com>.
Maintainer: Shachar Kaufman <shachark@post.tau.ac.il>
}
\references{
Heller, R., Heller, Y., & Gorfine, M. (2013). A consistent multivariate test of association based on ranks of distances. \emph{Biometrika}, 100(2), 503-510.

Heller, R., Heller, Y., Kaufman S., & Gorfine, M. (2014). Consistent distribution-free tests of association between univariate random variables. \emph{arXiv:1308.1559}.
}
%\seealso{
%}
\examples{
# See examples in the documentation for hhg.test and xdp.test.
}
