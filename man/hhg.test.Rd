\name{hhg.test}
\alias{hhg.test}

\title{Heller-Heller-Gorfine Tests of Independence}

\description{This function implements the Heller-Heller-Gorfine (HHG) test of independence between two random vectors (x and y) having arbitrary dimensions.}
\usage{
hhg.test(Dx, Dy, ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, 
         is.sequential = F, seq.total.nr.tests = 1, 
         seq.alpha.hyp = NULL, seq.alpha0 = NULL, seq.beta0 = NULL, seq.eps = NULL, 
         nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
}

\arguments{
  \item{Dx}{a symmetric matrix of doubles, where element [i, j] is a norm-based distance between the \code{i}'th and \code{j}'th \code{x} samples.}
  \item{Dy}{same as \code{Dx}, but for distances between \code{y}'s (the user may choose any norm when computing \code{Dx}, \code{Dy}).}
  \item{ties}{a boolean specifying whether ties in \code{Dx} and/or \code{Dy} exist and are to be properly handled (requires more computation).}
  \item{w.sum}{minimum expected frequency taken into account when computing the \code{sum.chisq} statistic (must be non-negative, contribution of tables having cells with smaller values will be truncated to zero).}
  \item{w.max}{minimum expected frequency taken into account when computing the \code{max.chisq} statistic (must be non-negative, contribution of tables having cells with smaller values will be truncated to zero).}
  \item{nr.perm}{number of permutations from which a p-value is to be estimated (must be non-negative). Can be specified as zero if only the observed statistics are wanted, without p-values. The actual number of permutations used may be slightly larger when using multiple processing cores. A Wald sequential probability ratio test is optionally implemented, which may push the p-value to 1 and stop permuting if it becomes clear that it is going to be high. See Details below.}
  \item{is.sequential}{boolean flag specifying whether Wald's sequential test is desired (see Details), otherwise a simple Monte-Carlo computation of \code{nr.perm} permutations is performed. When this argument is \code{TRUE}, either \code{seq.total.nr.tests} or (\code{seq.alpha.hyp}, \code{seq.alpha0}, \code{seq.beta0}, \code{seq.eps}) must be supplied by the user.}
  \item{seq.total.nr.tests}{the total number of hypotheses in the family of hypotheses simultaneously tested. When this optional argument is supplied, it is used to derive default values for the parameters of the Wald sequential test. The default derivation is done assuming a nominal \code{0.05} FDR level, and sets: \code{seq.alpha.hyp = 0.05 / max(1, log(seq.total.nr.tests))}, \code{seq.alpha0 = 0.05}, \code{seq.beta0 = min(0.01, 0.05 / seq.total.nr.tests)}, \code{seq.eps = 0.01}. Alternatively, one can specify their own values for these parameters using the following arguments.}
  \item{seq.alpha.hyp}{the nominal test size for this single test within the multiple testing procedure.}
  \item{seq.alpha0}{the nominal test size for testing the side null hypothesis of p-value > \code{seq.alpha.hyp}.}
  \item{seq.beta0}{one minus the power for testing the side null hypothesis of p-value > \code{seq.alpha.hyp}.}
  \item{seq.eps}{approximation margin around \code{seq.alpha.hyp} that defines the p-value regions for the side null \code{p > seq.alpha.hyp * (1 + seq.eps)} and side alternative \code{p < seq.alpha.hyp * (1 - seq.eps)}.}
  \item{nr.threads}{number of processing cores to use for p-value permutation. If left as zero, will try to use all available cores.}
  \item{tables.wanted}{boolean flag determining whether to output detailed local 2x2 contingency tables.}
  \item{perm.stats.wanted}{boolean flag determining whether to output statistics values computed for all permutations (representing null distributions).}
}

\details{
The HHG test (Heller et al., 2013) is a powerful nonparametric test for association (or, alternatively, independence) between two random vectors (say, x and y) of arbitrary dimensions. It is consistent against virtually any form of dependence, and has been shown to offer significantly more power than alternative approaches in the face of simple, and, more importantly, complex high-dimensional dependence. The test relies on norm-based distance metrics in x and (separately) in y. The choice of metric for either variable is up to the user (e.g. Euclidean, Mahalanobis, Manhattan, or whatever is appropriate for the data). The general implementation in hhg.test takes the distance matrices computed on an observed sample, Dx and Dy, and starts form there.

When enabled by \code{is.sequential}, Wald's sequential test is implemented as suggested by Fay et al. (2007) in order to reduce the \code{O(nr.perm * n^2 * log(n))} computational compelxity of the permutation test to a more managable size. Especially when faced with high multiplicity, say \code{M} simultaneous tests, the necessary number of iterations may be quite large. For example, if it is desired to control the family-wise error rate (FWER) at level \code{alpha} using the Bonferroni correction, one needs a p-value of \code{alpha / M} to establish significance. This seems to suggest that the minimum number of permutations required is \code{nr.perm = M / alpha}. However, if it becomes clear after a smaller number of permutations that the null cannot be rejected, no additional permutations are needed, and the p-value can be conservatively estimated as 1. Often, only a handful of hypotheses in a family are expected to be non-null. In this case the number of permutations for testing all hypotheses using Wald's procedure is expected to be much lower than the full \code{M^2 / alpha}. 

The target significance level of the sequential test is specified in the argument \code{seq.alpha.hyp}. It depends on the number of hypotheses \code{M}, and the type of multiplicity correction wanted. For the Bonferroni correction, the threshold is \code{alpha / M}. For the less conservative procedure of Benjamini & Hochberg (1995), it is \code{M1 * q / M}, where \code{q} is the desired false discovery rate (FDR), and \code{M1} is the (unknwon) number of true non-null hypotheses. Although \code{M1} is unknown, the investigator can sometimes estimate it conservatively (e.g., if at most 0.02 of the hypotheses are expected to be non-null, set \code{M1 = 0.02 * M}).
}

\value{
Four statistics described in the original HHG paper are returned:

\code{sum.chisq} - sum of Pearson chi-squared statistics from the 2x2 contingency tables considered
\code{sum.lr} - sum of liklihood ratio ("G statistic") values from the 2x2 tables
\code{max.chisq} - maximum Pearson chi-squared statistic from any of the 2x2 tables
\code{max.lr} - maximum G statistic from any of the 2x2 tables

If \code{nr.perm > 0}, then estimated permutation p-values are returned as:

\code{perm.pval.hhg.sc}
\code{perm.pval.hhg.sl}
\code{perm.pval.hhg.mc}
\code{perm.pval.hhg.ml}

In order to give information that may help localize where in the support of the distributions of \code{x} and \code{y} there is departure from independence, if \code{tables.wanted} is true, the 2x2 tables themselves are provided in:

\code{extras.hhg.tbls}

This is a \code{n^2} by 4 matrix, whose columns are \code{A11}, \code{A12}, \code{A21}, \code{A22} as denoted in the original HHG paper. Row r of the matrix corresponds to S_ij in the same paper, where \code{i = 1 + floor((r - 1) / n)}, and \code{j = 1 + ((r - 1) \%\% n)}. Since \code{S_ij} is never computed for \code{i == j}, rows \code{(0:(n - 1)) * n + (1:n)} contain \code{NA}s on purpose.

Finally, as a means of estimating the null distributions of computed statistics, if \code{perm.stats.wanted} is true, the statistics computed for every permutation of the data performed during testing is outputted as:

\code{extras.perm.stats}

A data.frame with one variable per statistic and one sample per permutation.
}

\references{
Heller, R., Heller, Y., & Gorfine, M. (2013). A consistent multivariate test of association based on ranks of distances. \emph{Biometrika}, 100(2), 503-510.

Fay, M., Kim., H., & Hachey, M. (2007). On using truncated sequential probability ratio test boundaries for Monte Carlo implementation of hypothesis tests. \emph{Journal of Computational and Graphical Statistics}, 16(4), 946-967.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series B}, 57, 289-300.

Hall, P. & Tajvidi, N. (2002). Permutation tests for equality of distributions in high-dimensional settings. \emph{Biometrika}, 89(2), 359-374.

Szekely, G., Rizzo, M., & Bakirov, N. (2007). Measuring and testing independence by correlation of distances. \emph{The Annals of Statistics}, 35, 2769-2794.
}

\author{
  Shachar Kaufman, based in part on an earlier version by Ruth Heller and Yair Heller.
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
