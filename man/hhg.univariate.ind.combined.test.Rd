\name{hhg.univariate.ind.combined.test}
\alias{hhg.univariate.ind.combined.test}


\title{The combined independence test statistics for all partition sizes}

\description{A Class of Consistent Univariate Distribution-Free Tests of Independence - function for computation of combined test}

\usage{
hhg.univariate.ind.combined.test(X,Y=NULL,NullTable=NULL,mmin=2,
  mmax=max(floor(sqrt(length(X))/2),2),variant='ADP',aggregation.type='sum',
  score.type='LikelihoodRatio', w.sum = 0, w.max = 2 ,combining.type='MinP',nr.perm=100)
}

\arguments{
  \item{X}{a numeric vector with observed \code{X} values,  or the test statistic as output from \code{\link{hhg.univariate.ind.stat}}.}
  \item{Y}{a numeric vector with observed \code{Y} values.  Leave as Null if the input to \code{X} is the test statistic. }
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website or computed by the function \code{\link{hhg.univariate.ind.nulltable}}.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is half the square root of the number of observations. For a max aggregation.type, this parameter cannot be more than 2 for the ADP varaint and 4 for DDP variant}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. }
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
  \item{combining.type}{a character string specifying the combining type, must be one of \code{"MinP"} (default), \code{"Fisher"}, or \code{"both"}.}
  \item{nr.perm}{The number of permutations for the null distribution. Ignored if \code{NullTable} is non-null.}
}

\details{
  The test statistic and p-value of the recommended independence test between two univariate random variables in Heller et al. (2014). The default combining type in the minimum p-value, so the test statistic is the minimum p-value over the range of partition sizes m from \code{mmin} to \code{mmax}, where the p-value for a fixed partition size m is defined by the aggregation type and score type . 
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  \item{MinP}{The test statistic when the combining type is \code{"MinP"}.}
  
  \item{MinP.pvalue}{The p-value when the combining type is \code{"MinP"}.}
  
  \item{MinP.m.chosen}{The partition size m for which the p-value was the smallest.}
  
  \item{Fisher}{The test statistic when the combining type is \code{"Fisher"}.}
  
  \item{Fisher.pvalue}{The p-value when the combining type is \code{"Fisher"}.}
  
  \item{m.stats}{The statistic for each m in the range \code{mmin} to \code{mmax}.}
  
  \item{pvalues.of.single.m}{The p-values for each m in the range \code{mmin} to \code{mmax}.}
  
  \item{generated_null_table}{The null table object. Null if \code{NullTable} is non-null.}
  
  \item{stat.type}{"Independence-Combined"}
  
  \item{variant}{a character string specifying the partition type used in the test, one of \code{"ADP"} or \code{"DDP"}.}
   
  \item{aggregation.type}{a character string specifying the aggregation type used in the , one of \code{"sum"} or \code{"max"}.}
  
  \item{score.type}{a character string specifying the score typeused in the test, one of \code{"LikelihoodRatio"} or \code{"Pearson"}. }
  
  \item{mmax}{The maximum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{mmin}{The minimum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}.}
  

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

#I) Perform MinP & Fisher Tests - without existing null tables.
#Null tables are generated by the test function.
#using partitions sizes up to 5
results = hhg.univariate.ind.combined.test(X,Y,nr.perm = 100,mmax=5)
results


#The null table can then be accessed.
generated.null.table = results$generated_null_table


#II) Perform MinP & Fisher Tests - with existing null tables. 

#create null table for aggregation by summation (on ADP), with partitions sizes up to 5: 
ADP.null = hhg.univariate.ind.nulltable(N,mmax=5)

#create a null table, using aggregation by summation over DDP partitions,
#with partitions sizes up to 5, over Pearson scores,
#with 1000 bootstrap repetitions.
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
score.type = 'Pearson', nr.replicates = 1000)

MinP.ADP.existing.null.table = hhg.univariate.ind.combined.test(X,Y, NullTable = ADP.null)
#Results 
MinP.ADP.existing.null.table

#using the other null table (DDP variant, with pearson scores):
MinP.DDP.existing.null.table = hhg.univariate.ind.combined.test(X,Y, NullTable = DDP.null)

MinP.DDP.existing.null.table

# combined test can also be performed by using the test statistic.
ADP.statistic = hhg.univariate.ind.stat(X,Y,mmax=5)
MinP.using.statistic.result = hhg.univariate.ind.combined.test(ADP.statistic, NullTable = ADP.null)
# same result as above (as MinP.ADP.result.using.existing.null.table$MinP.pvalue)
MinP.using.statistic.result$MinP.pvalue

}


}
