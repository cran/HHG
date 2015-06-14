\name{hhg.univariate.ind.nulltable}
\alias{hhg.univariate.ind.nulltable}


\title{The independence test null tables for all partition sizes}

\description{Functions for creating null table objects, for the omnibus distribution-free test of independence between two univariate random variables, as described in Heller et al. (2014). To be used for the p-value computation, see examples in \code{\link{hhg.univariate.ind.pvalue}}.}

\usage{
hhg.univariate.ind.nulltable(size,mmin=2,mmax=max(floor(sqrt(size)/2),2),
variant='ADP',aggregation.type='sum',score.type='LikelihoodRatio',
w.sum = 0, w.max = 2,nr.replicates=1000,keep.simulation.data=F)
}

\arguments{
  \item{size}{The sample size}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is half the square root of the number of observations.}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default) or \code{"max"}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default) or \code{"Pearson"}. }
   \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
  \item{nr.replicates}{The number of permutations for the null distribution.}
  \item{keep.simulation.data}{TRUE/FALSE. If TRUE, then in addition to the sorted statistics per column, the original matrix of size nr.replicates by mmax-mmin+1 is also stored.}
}

\details{
 In order to compute the null distributions for a test statistic  (with a specific aggregation and score type, and all partition sizes), the only necessary information is the sample size. The accuracy of the quantiles of the null distribution depend on the number of replicates used for constructing the null tables. The necessary accuracy depends on the threshold used for rejection of the null hypotheses.
}

\value{
   \item{m.stats}{If keep.simulation.data= TRUE, \code{m.stats} a matrix with \code{nr.replicates} rows and  \code{mmax-mmin+1} columns of null test statistics.}
  
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
#Testing for independance, sample size = 35
N=35

#null table for aggregation by summation (on ADP): 
ADP.null = hhg.univariate.ind.nulltable(N)

#create a null table, using aggregation by summation over DDP partitions,
#with partitions sizes up to 5, over Pearson scores, with 1000 bootstrap repetitions.
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
score.type = 'Pearson', nr.replicates = 1000)
}


}
