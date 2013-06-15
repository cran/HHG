\name{HHG}
\alias{hhg.test}

\title{Heller-Heller-Gorfine Tests of Independence}

\description{These functions perform Heller-Heller-Gorfine (HHG) tests. Implemented are tests of independence between two random vectors (x and y).}
\usage{
hhg.test(Dx, Dy, ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, 
         total.nr.tests = 1, is.sequential = T, 
         alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
         nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
}

\arguments{
  \item{Dx}{a symmetric matrix of doubles, where element [i, j] is a norm-based distance between x_i and x_j.}
  \item{Dy}{same as \code{Dx}, but for distances between y_i and y_j (user may choose any norm when computing \code{Dx}, \code{Dy}).}
  \item{ties}{a boolean specifying whether ties in Dx and/or Dy exist and are to be properly handled (requires more computation).}
  \item{w.sum}{minimum expected frequency taken into account when computing the sum_chi statistic (must be non-negative, contribution of cells with smaller values will be truncated to zero).}
  \item{w.max}{minimum expected frequency taken into account when computing the max_chi statistic (must be non-negative, contribution of cells with smaller values will be truncated to zero).}
  \item{nr.perm}{number of permutations from which a p-value is to be estimated (must be non-negative). Can be specified as zero if only the observed statistics are wanted, without p-values. The actual maximum number used may be slightly larger when using multiple processing cores. A Wald sequential probability ratio test is optionally implemented (ans is the default) which may push the p-value to 1 if it becomes clear it is going to be high. See details below.}
  \item{is.sequential}{boolean flag whether Wald's sequential test is desired, otherwise a simple Monte-Carlo computation of nr.perm permutations is performed.}
  \item{total.nr.tests}{the total number of hypotheses in the family of hypotheses simultaneously tested. This is used to derive default values for the parameters of the Wald sequential test, assuming a nominal 0.05 FDR level: \code{alpha.hyp = 0.05 / max(1, log(total.nr.tests)), alpha0 = 0.05, beta0 = min(0.01, 0.05 / total.nr.tests), eps = 0.01} (in which case, the user does not have to supply them). This is only relevant if \code{is.sequential} is true.}
  \item{alpha.hyp}{the nominal test size for this single test within the multiple testing procedure.}
  \item{alpha0}{the nominal test size for testing the side null hypothesis of p-value > \code{alpha.hyp}.}
  \item{beta0}{one minus the power for testing the side null hypothesis of p-value > \code{alpha.hyp}.}
  \item{eps}{approximation margin around alpha.hyp that defines the p-value regions for the side null \code{p > alpha.hyp * (1 + eps)} and side alternative \code{p < alpha.hyp * (1 - eps)}.}
  \item{nr.threads}{number of processing cores to use for p-value permutation. If left as zero, will try to use all available cores.}
  \item{tables.wanted}{boolean flag determining whether to output detailed local 2x2 contingency tables.}
  \item{perm.stats.wanted}{boolean flag determining whether to output statistics values computed for all permutations (representing null distributions).}
}

\details{
The HHG test (Heller et al. 2012) is a powerful nonparametric test for association (or, alternatively, independence) between two random vectors (say, x and y) of arbitrary dimensions. It is consistent against any form of dependence, and has been shown to offer significantly more power than alternative approaches in the face of simple, and more importantly, complex high-dimensional dependence. The test relies on norm-based distance metrics in x and (separately) in y. The choice of metric for either variable is up to the user (e.g. Euclidean, Mahalanobis, Manhattan, or whatever is appropriate for the data). The general implementation in hhg.test takes the distance matrices computed on an observed sample, Dx and Dy, and starts form there.

When enabled by \code{is.sequential}, Wald's sequential test is implemented as suggested by Fay et al. (2007) in order to reduce the O(nr.perm * n^2 * log(n)) computational compelxity of the permutation test to a more managable size. When faced with multiplicity, the necessary number of iterations may be quite large. For example, if it is desired to control the FWER at level 0.05 with a Bonferroni correction for multiplicity, the minimum number of iterations to establish significance of a p-value is monte = 20 * total.nr.tests. However, if the hypothesis is deemed uninteresting after a smaller number of iteration, no additional permutations are needed, and the p-value can be conservatively estimated as 1. Often, only a handful of hypotheses are expected to be non-null. In this case the complexity for all hypotheses is expected to be much lower than 20 * (total.nr.tests)^2. The threshold for decideding which hypotheses are interesting depends on the total number of hypotheses simultaneously examined. For a Bonferroni correction, this threshold is alpha/total.nr.tests, where alpha is the desired FWER level (say alpha = 0.05). For the less conservative procedure of Benjamini & Hochberg (1995), this threshold is M1 * q / total.nr.tests, where q is the desired FDR level (say q = 0.05), and M1 is the unknwon number of true non-null hypothesese (i.e., dependent associations). Although M1 is unknown, it may conservatively be estimated by the investigator (e.g., if at most 0.02 of the hypotheses are expected to be non-null, set M1 = 0.02 * total.nr.tests).  
}

\value{
Four statistics described in the original HHG paper are returned:

\code{sum.chisq} - sum of Pearson chi-squared statistics from the 2x2 contingency tables considered
\code{sum.lr} - sum of liklihood ratio ("G statistic") values from the 2x2 tables
\code{max.chisq} - maximum Pearson chi-squared statistic from any of the 2x2 tables
\code{max.lr} - maximum G statistic from any of the 2x2 tables

Their estimated p-values are returned as:

\code{perm.pval.hhg.sc}
\code{perm.pval.hhg.sl}
\code{perm.pval.hhg.mc}
\code{perm.pval.hhg.ml}

In order to give information that may help localize where in the support of the distributions of \code{x} and \code{y} there is departure from independence, if \code{tables.wanted} is true, the 2x2 tables themselves are provided in:

\code{extras.hhg.tbls}

This is a n^2 by 4 matrix, whose columns are A11, A12, A21, A22 as denoted in the original HHG paper. Row r of the matrix corresponds to S_ij in the same paper, where \code{i = 1 + floor((r - 1) / n)}, and \code{j = 1 + ((r - 1) \%\% n)}. Since S_ij is never computed for \code{i == j}, rows \code{(0:(n - 1)) * n + (1:n)} contain \code{NA}s on purpose. The only other case where NAs will occur are for the 2 and K-sample tests, where only one table is given for any x-tied samples (the other tables at indices with the same x value are redundant).

Finally, as a means of estimating the null distributions of computed statistics, if \code{perm.stats.wanted} is true, the statistics computed for every permutation of the data performed during testing is outputted as:

\code{extras.perm.stats}

A data.frame with one variable per statistic and one sample per permutation.
}

\references{
Heller R., Heller Y., and Gorfine M. (2012). A consistent multivariate test of association based on ranks of distances. \emph{arXiv:1201.3522v1}.

Fay, M., and Kim., H., and Hachey, M. (2007). On Using Truncated Sequential Probability Ratio Test Boundaries for Monte Carlo Implementation of Hypothesis Tests. \emph{Journal of Computational and Graphical Statistics}, 16(4): 946-967.

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series} B, \bold{57}, 289-300.
}

\author{
  Shachar Kaufman, based on an earlier version by Ruth Heller and Yair Heller.
}

\note{
  The computational complexity of the test is \code{n^2*log(n)}, where \code{n} is the number of samples. Thus, when the sample size is large, computing the test for many permutations may take a long time.
}

\examples{
## 1. The test of independence

## 1.1. A non-null univariate example

## Generate some data from the Circle example
n = 50
X = hhg.example.datagen(n, 'Circle')
plot(X[1,], X[2,])

## Compute distance matrices, on which the HHG test will be based
Dx = as.matrix(dist((X[1,]), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((X[2,]), diag = TRUE, upper = TRUE))

## Compute HHG statistics, and p-values using 1000 random permutations
hhg = hhg.test(Dx, Dy, nr.perm = 1000)

## Print the sum-chisq statistic and its permutation p-value
hhg$sum.chisq
hhg$perm.pval.hhg.sc

## 1.2. A null univariate example

n = 50
X = hhg.example.datagen(n, '4indclouds') 

Dx = as.matrix(dist((X[1,]), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((X[2,]), diag = TRUE, upper = TRUE))

hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg$sum.chisq
hhg$perm.pval.hhg.sc

hhg$sum.lr
hhg$perm.pval.hhg.sl

## 1.3. A multivariate example
library(MASS)

n = 50
p = 5
x = t(mvrnorm(n, rep(0, p), diag(1, p)))
y = log(x ^ 2)
Dx = as.matrix(dist((t(x)), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((t(y)), diag = TRUE, upper = TRUE))

hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg$sum.chisq
hhg$perm.pval.hhg.sc
}
