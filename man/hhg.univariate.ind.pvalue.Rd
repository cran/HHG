\name{hhg.univariate.ind.pvalue}
\alias{hhg.univariate.ind.pvalue}


\title{The p-value computation for the test of independence using a fixed partition size}

\description{The p-value computation for the distribution free test of independence between two univariate random variables of Heller et al. (2014) ,using a fixed partition size \code{m}.}

\usage{
hhg.univariate.ind.pvalue(statistic, NullTable, m=min(statistic$mmax,4))
}

\arguments{
  \item{statistic}{The value of the computed statistic  by the function \code{\link{hhg.univariate.ind.stat}}. The statistic object includes the score type (one of \code{"LikelihoodRatio"} or \code{"Pearson"}), and the aggregation type (one of \code{"sum"} or \code{"max"}).}
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website (http://www.math.tau.ac.il/~ruheller/Software.html) or computed by the function
  
  \code{\link{hhg.univariate.ind.nulltable}}. See \code{vignette('HHG')} for a method of computing null tables on multiple cores.}
  \item{m}{The partition size.}
}



\details{
  For the test statistic, the function extracts the fraction of observations in the null table that are at least as large as the test statistic, i.e. the p-value. 
}

\value{
  The p-value. 

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


#I) Computing test statistics , with default parameters:

#statistic:
hhg.univariate.ADP.Likelihood.result = hhg.univariate.ind.stat(X,Y)
hhg.univariate.ADP.Likelihood.result

#null table:
ADP.null = hhg.univariate.ind.nulltable(N)
#pvalue:
hhg.univariate.ind.pvalue(hhg.univariate.ADP.Likelihood.result, ADP.null)

#II) Computing test statistics , with summation over Data Derived Partitions (DDP),
#using Pearson scores, and partition sizes up to 5:

#statistic:
hhg.univariate.DDP.Pearson.result = hhg.univariate.ind.stat(X,Y,variant = 'DDP',
  score.type = 'Pearson', mmax = 5)
hhg.univariate.DDP.Pearson.result

#null table:
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
  score.type = 'Pearson', nr.replicates = 1000)
  
#pvalue , for different partition size:
hhg.univariate.ind.pvalue(hhg.univariate.DDP.Pearson.result, DDP.null, m =2)
hhg.univariate.ind.pvalue(hhg.univariate.DDP.Pearson.result, DDP.null, m =5)

}



}
