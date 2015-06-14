\name{hhg.univariate.ks.nulltable}
\alias{hhg.univariate.ks.nulltable}


\title{The K-sample test null tables for all partition sizes}

\description{Functions for creating null table objects, for the omnibus distribution-free test of equality of distributions among K groups, as described in Heller et al. (2014). To be used for the p-value computation, see examples in \code{\link{hhg.univariate.ks.pvalue}}. }

\usage{
hhg.univariate.ks.nulltable(group.sizes,mmin=2,
mmax=max(4,round(min(group.sizes)/3)), aggregation.type='sum',
score.type='LikelihoodRatio',nr.replicates=1000,keep.simulation.data=F)
}

\arguments{
  \item{group.sizes}{the number of observations in each group.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is 1/3 the number of observations in the smallest group.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default) or \code{"max"}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default) or \code{"Pearson"}.}
  \item{nr.replicates}{The number of permutations for the null distribution.}
  \item{keep.simulation.data}{TRUE/FALSE. If TRUE, then in addition to the sorted statistics per column, the original matrix of size nr.replicates by mmax-mmin+1 is also stored.}
}

\details{
  In order to compute the null distributions for a test statistic  (with a specific aggregation and score type, and all partition sizes), the only necessary information is the group sizes. The accuracy of the quantiles of the null distribution depend on the number of replicates used for constructing the null tables. The necessary accuracy depends on the threshold used for rejection of the null hypotheses.  
}

\value{
  \item{m.stats}{If keep.simulation.data= TRUE, a matrix with \code{nr.replicates} rows and  \code{mmax-mmin+1} columns of null test statistics.}
  
  \item{univariate.object}{A useful format of the null tables for computing p-values efficiently.}

}

\references{
Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2014). Consistent distribution-free K-sample and independence tests for univariate random variables \emph{arXiv:1410.6758}.
}

\author{
  Barak Brill and Shachar Kaufman.
}



\examples{
\dontrun{
#Testing for the difference between two groups, each from a normal mixture:
N0=30
N1=30

#null table for aggregation by summation: 
sum.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), nr.replicates=100)
#default nr. of replicates is 1000,
#but may take several seconds. For illustration only, we use 100 replicates,
#but it is highly recommended to use at least 1000 in practice. 

#null table for aggregation by maximization: 
max.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), aggregation.type = 'max',
  score.type='LikelihoodRatio', mmin = 3, mmax = 5, nr.replicates = 100)
#default nr. of replicates is 1000, but may take several seconds. For illustration only,
#we use 100 replicates, but it is highly recommended to use at least 1000 in practice.

}

}
