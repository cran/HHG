\name{hhg.example.datagen}
\alias{hhg.example.datagen}
\title{
A set of example data generators used to demonstrate the HHG test.
}
\description{
Six examples (Circle, Diamond, Parabola, 2Parabolas, W, 4indclouds) are taken from Newton's introduction to the discussion of the Energy Test in The Annals of Applied Statistics (2009). These are simple univariate dependence structures (or independence, in the latter case) used to demonstrate the tests of independece. The remaining examples (TwoClassUniv, FourClassUniv, TwoClassMultiv) generate data suitable for demonstrating the k-sample test (and in particular, the two-sample test).
}
\usage{
hhg.example.datagen(n, example)
}
\arguments{
  \item{n}{The desired sample size}
  \item{example}{The choice of example}
}
\value{
  For \code{example} in \{Circle, Diamond, Parabola, 2Parabolas, W, and 4indclouds\}, a matrix of two rows is returned, one row per variable. Columns are i.i.d. samples. Given these data, we would like to test whether the two variables are statistically independent. Except for the 4indclouds case, all examples in fact have variables that are dependent.
  When \code{example} is one of \{TwoClassUniv, FourClassUniv, TwoClassMultiv\}, a list is returned with elements \code{x} and \code{y}. \code{y} is a vector with values either 0 or 1 (for TwoClassUniv and TwoClassMultiv) or in 0:3 for (for FourClassUniv). \code{x} is a real valued random variable (TwoClassUniv and FourClassUniv) or vector (TwoClassMultiv) which is not independent of \code{y}.
}
\references{
Newton, M.A. (2009). Introducing the discussion paper by Szekely and Rizzo. \emph{Annals of applied statistics}, 3 (4), 1233-1235.
}
\author{
Shachar Kaufman and Ruth Heller
}
\examples{
X = hhg.example.datagen(50, 'Diamond')
plot(X[1,], X[2,])

X = hhg.example.datagen(50, 'FourClassUniv')
plot(X)
}