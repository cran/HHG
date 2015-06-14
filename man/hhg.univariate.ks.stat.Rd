\name{hhg.univariate.ks.stat}
\alias{hhg.univariate.ks.stat}


\title{The K-sample test statistics for all partition sizes}

\description{These statistics are used in the omnibus distribution-free test of equality of distributions among K groups, as described in Heller et al. (2014). }

\usage{
hhg.univariate.ks.stat(x, y,aggregation.type='sum',score.type='LikelihoodRatio',
mmax = max(4,round(min(table(y))/3)),mmin=2)
}

\arguments{
  \item{x}{a numeric vector of data values. Tied observations are broken at random.}
  \item{y}{for \code{k} groups, a vector of integers with values \code{0:(k-1)} which specify the group each observation belongs to.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}} 
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. }
  \item{mmax}{The maximum partition size of the ranked observations, default value is 1/3 the number of observations in the smallest group.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
}

\details{
  For each partition size \eqn{m= mmin,\ldots,mmax}, the function computes the scores in each of the paritions (according to score type), and aggregates all scores according to the aggregation type (see details in Heller et al. , 2014). If the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}, then the computed statistic will be in \code{statistic}, otherwise the computed statistics will be in the appropriate subset of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}.
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  
  \item{statistic}{The value of the computed statistic if the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}. One of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}. }
  
  \item{sum.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered, divided by the total number of observations.} 
  
  \item{sum.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all LikelihoodRatio statistics from all the \eqn{K X m} contingency tables considered, divided by the total number of observations.} 
  
 \item{max.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered.} 
  
  \item{max.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered.} 
  
  \item{type}{"KSample".}
  
  \item{stat.type}{"KSample".}
  
  \item{size}{A vector of size K of the ordered group sample sizes.}
  
  \item{score.type}{The input \code{score.type}.}
  
  \item{aggregation.type}{The input \code{aggregation.type}.}
  
  \item{mmin}{The input \code{mmin}.}
  
  \item{mmax}{The input \code{mmax}.}

}

\references{
Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2014). Consistent distribution-free K-sample and independence tests for univariate random variables \emph{arXiv:1410.6758}.
}

\author{
  Barak Brill and Shachar Kaufman.
}



\examples{
#Example of computing the test statisics for data from a two-sample problem:

#Two groups, each from a different normal mixture:
X = c(c(rnorm(25,-2,0.7),rnorm(25,2,0.7)),c(rnorm(25,-1.5,0.5),rnorm(25,1.5,0.5)))
Y = (c(rep(0,50),rep(1,50)))
plot(Y,X)


#I) Computing test statistics , with default parameters:
hhg.univariate.Sm.Likelihood.result = hhg.univariate.ks.stat(X,Y)

hhg.univariate.Sm.Likelihood.result

#II) Computing test statistics , with max aggregation type:
hhg.univariate.Mm.likelihood.result = hhg.univariate.ks.stat(X,Y,aggregation.type = 'max')

hhg.univariate.Mm.likelihood.result


}
