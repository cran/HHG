\name{XDP}
\alias{xdp.test}

\title{A Class of Consistent Univariate Distribution-Free Tests of Independence}

\description{This function computes test statistics for the tests of independence between two random variables (\code{X} and \code{Y}) presented in Kaufman et al. (2013). Several specialized variants of data-driven partitions (DDP) and all data partitions (ADP) tests are implemented.}

\usage{
xdp.test(x, y, variant = 'DDP', K = 3, correct.mi.bias = F, w.sum = 0, w.max = 2)
}

\arguments{
  \item{x}{a numeric or ordered factor vector with observed \code{X} values.}
  \item{y}{a numeric or ordered factor vector with observed \code{Y} values, of the same length as \code{x}.}
  \item{variant}{either 'DDP', 'ADP'.}
  \item{K}{the dimension of the KxK partition. Must be between 2 and the square root of the length of the input vectors.}
  \item{correct.mi.bias}{a boolean specifying whether to apply the Miller-Madow bias correction that can be useful for mutual information estimation.}
  \item{w.sum}{minimum expected frequency taken into account when computing the sum_chi statistic (must be non-negative, contribution of tables having cells with smaller values will be truncated to zero). Note this only effects the specialized algorithms for small \code{K} as described below.}
  \item{w.max}{minimum expected frequency taken into account when computing the max_chi statistic (must be non-negative, contribution of tables having cells with smaller values will be truncated to zero). Note this only effects the specialized algorithms for small \code{K} as described below.}
}

\details{
This is an implementation of the DDP (data-driven partition) / ADP (all data partitions) class of univariate distribution-free tests (Kaufman et al., 2013). This class extends Hoeffding's test by considering data-driven (or all) partitions of an arbitrary order \code{KxK} of the paired sample (or rank) sapce. Being distribution-free, it does not require resampling for computing p-values, making it considerably faster to apply than the HHG test (see example below).

Specialized algorithms are used for DDP with \code{K=2,3,4}, and ADP with \code{K=2}, which are of \code{O(n^1), O(n^2), O(n^3)} and \code{O(n^2)} time complexity, respectively (with \code{n = length(x)}). These algorithms compute, in addition to the sum-aggregated statistics which are available for any \code{K}, statistics based on max aggregation. These algorithms take \code{w.sum} and \code{w.max} into account.

For generating p-values, either simulate your own null tables as shown in the example below, or, if your sample size (\code{length(x)}) is available, use the pre-computed tables linked from Dr. Ruth Heller's website at \url{http://www.math.tau.ac.il/~ruheller/Software.html}. See the README file bundled with the tables for current details on availability and storage format. At the time of writing tables are available for sample sizes 30, 40, 50, 100, and 300. For each sample size, DDP and ADP tests with \code{K = 2:floor(sqrt(length(x)))} were computed. The tables contain 1e6 replicates (i.e., support p-values down to 1e-6), except for sample size 300 where currently we have 2e4 samples.
}

\value{
The specialized cases return four statistics:

\code{sum.chisq} - normalized sum of Pearson chi-squared statistics from the KxK contingency tables considered.

\code{sum.lr} - normalized sum of liklihood ratio ("G statistic") values from the KxK tables.

\code{max.chisq} - maximum Pearson chi-squared statistic from any of the KxK tables.

\code{max.lr} - maximum G statistic from any of the KxK tables.

Sum statistics are normalized by the total number of partitions multiplied by \code{length(x)}, which is the scale for measuring mutual information (in natural digits, aka 'nats').

As described in Kaufman et al. (2013), DDP with \code{K > 4} and ADP with \code{K > 2} resort to by-cell rather than by-partition enumeration, and thus cannot currently return 'max' statistics (only 'sum' statistics are returned).
}

\references{
Heller, R., Heller, Y., Kaufman S., & Gorfine, M. (2014). Consistent distribution-free tests of association between univariate random variables. \emph{arXiv:1308.1559}.
}

\author{
  Shachar Kaufman.
}

\note{
The tests are consistent for continuous random variables, i.e. when there are no ties. Any ties that do exist in \code{x} and \code{y} are automatically resolved by random assignment of ranks (see the built-in function \code{rank}).
}

\examples{

# Background: Similar expression (activity) patterns of genes indicate that 
# those genes share a common function. It is thus possible to infer the function
# of genes where it is currently unknown, or to learn about pathways involving 
# multiple genes working in concert. Hopefully this knowledge can then be used 
# to develop better crops, new treatments to disease, and so on.

# In this example we use the DDP test to find pairs of genes with associated
# expression patterns in the data of Hughes et al. (2000). Specifically, we will
# use the 3x3 DDP test variant based on summation of the likelihood ratio.

data(hughes)
nr.genes = nrow(hughes)
nr.samples = ncol(hughes)

# We want to test each pair of genes for association.
nr.tests = choose(nr.genes, 2)

# Estimate the critical test statistic value by Monte-Carlo simulation of the
# null distribution. We use enough null replicates in order to ensure accurate
# estimation.
bonferroni.level = 0.05 / nr.tests
nr.null.replicates = 100 / bonferroni.level

if (0) {
  # This takes a while to compute but only needs to be done once (and is easy to
  # parallelize over many CPUs).
  ddp.nullsim = rep(Inf, nr.null.replicates)
  x = 1:nr.samples # it is enough to permute the y vector
  
  for (i in 2:nr.null.replicates) {
    y = sample(nr.samples)
    ddp.nullsim[i] = xdp.test(x, y, variant = 'DDP', K = 3)$sum.lr
  }
  
  ddp.critical.value = quantile(ddp.nullsim, probs = 1 - bonferroni.level)
} else {
  # In the interest of saving time, this is what we get
  # NOTE: Big Monte-Carlo simulated null distributions are available (see documentation)
  ddp.critical.value = 0.01569162
}

# Now we compute the observed test statistics for all pairs
compute.xdp.for.pair = function(pair.idxs) {
  x = as.numeric(hughes[pair.idxs[1], ])
  y = as.numeric(hughes[pair.idxs[2], ])
  return (xdp.test(x, y, variant = 'DDP', K = 3)$sum.lr)
}

ddp.obs = combn(nr.genes, 2, compute.xdp.for.pair, simplify = TRUE)

# And identify any significantly associated pairs (at the conservative 
# Bonferroni level maintaining the FWER at 0.05)
ddp.significant = (ddp.obs > ddp.critical.value)

summary(ddp.significant)

# Result: There are 26 pairs out of the 45 performed here, where there is 
# evidence of associated patterns! Had we analyzed this data with Pearson or 
# Spearman correlation based testing, we would find only a small fraction of 
# these associated pairs (see Kaufman et al., 2013).
}
