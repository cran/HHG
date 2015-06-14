\name{hhg.univariate.ind.stat}
\alias{hhg.univariate.ind.stat}


\title{The independence test statistics for all partition sizes}

\description{These statistics are used in the omnibus distribution-free test of independence between two univariate random variables, as described in Heller et al. (2014).}

\usage{
hhg.univariate.ind.stat(x, y, variant = 'ADP',aggregation.type='sum'
  ,score.type='LikelihoodRatio',mmax = max(floor(sqrt(length(x))/2),2),
  mmin =2 , w.sum = 0, w.max = 2)
}

\arguments{
  \item{x}{a numeric vector with observed \code{X} values (tied observations are broken at random).}
  \item{y}{a numeric vector with observed \code{Y} values (tied observations are broken at random).}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. }
  \item{mmax}{The partition size of the ranked observations. The default size is half the square root of the number of observations}
  \item{mmin}{The partition size of the ranked observations. The default size is half the square root of the number of observations}
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
}

\details{
  For each partition size \eqn{m= mmin,\ldots,mmax}, the function computes the scores in each of the paritions (according to score type), and aggregates all scores according to the aggregation type (see details in Heller et al. , 2014). If the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}, then the computed statistic will be in \code{statistic}, otherwise the computed statistics will be in the appropriate subset of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}. Note that if the variant is \code{"ADP"}, all partition sizes are computed together in O(N^4), so the score computational complexity is O(N^4). For \code{"DDP"} and mmax>4,the score computational complexity is O(N^4)*(mmax-mmin+1).  
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  
  \item{statistic}{The value of the computed statistic if the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}. One of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}.}
  
  \item{sum.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered, divided by the total number of observations.}
  
  \item{sum.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all LikelihoodRatio statistics from all the \eqn{m X m} contingency tables considered, divided by the total number of observations.}
  
  \item{max.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered.}
  
  \item{max.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered.}
  
  \item{type}{"Independence"}
  
  \item{stat.type}{"Independence-Stat"}
  
  \item{size}{The sample size}
  
  \item{score.type}{The input \code{score.type}.}
  
  \item{aggregation.type}{The input \code{aggregation.type}.}
  
  \item{mmin}{The input \code{mmin}.}
  
  \item{mmax}{The input \code{mmax}.}
  
  \item{additional}{A vector with the input \code{w.sum} and \code{w.max}.}
 

}

\references{
Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2014). Consistent distribution-free K-sample and independence tests for univariate random variables \emph{arXiv:1410.6758}.
}

\author{
  Barak Brill and Shachar Kaufman.
}


\examples{

\dontrun{
N = 35
data = hhg.example.datagen(N, 'Parabola')
X = data[1,]
Y = data[2,]
plot(X,Y)


#I) Computing test statistics , with default parameters(ADP statistic):

hhg.univariate.ADP.Likelihood.result = hhg.univariate.ind.stat(X,Y)

hhg.univariate.ADP.Likelihood.result

#II) Computing test statistics , with summation over Data Derived Partitions (DDP),
#using Pearson scores, and partition sizes up to 5:

hhg.univariate.DDP.Pearson.result = hhg.univariate.ind.stat(X,Y,variant = 'DDP',
  score.type = 'Pearson', mmax = 5)
hhg.univariate.DDP.Pearson.result
}

}
