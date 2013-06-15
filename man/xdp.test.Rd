\name{XDP}
\alias{xdp.test}

\title{A Class of Consistent Univariate Distribution-Free Tests of Independence}

\description{This function computes test statistics for the tests of independence between two random variables (\code{X} and \code{Y}) presented in Kaufman et al. (2013). Several specialized variants of data-driven partitions (DDP) and all data partitions (ADP) tests are implemented.}

\usage{
xdp.test(x, y, variant = 'ppr.33.obs', K = 3, correct.mi.bias = F)
}

\arguments{
  \item{x}{a numeric or ordered factor vector with observed \code{X} values.}
  \item{y}{a numeric or ordered factor vector with observed \code{Y} values, of the same length as \code{x}.}
  \item{variant}{either 'spr.obs', 'spr.all', 'ppr.22.obs', 'ppr.33.obs', 'ppr.22.all', 'ppr.33.all', 'tpr.obs', 'tpr.all', 'ddp.obs', or 'ddp.all'. See details below.}
  \item{K}{the dimension of the KxK partition. Relevant only for the 'ddp' variants, and must be between 2 and the square root of the length of the input vectors.}
  \item{correct.mi.bias}{a boolean specifying whether to apply the Miller-Madow bias correction that can be useful for mutual information estimation.}
}

\details{
This is an implementation of the DDP/ADP class of univariate distribution-free tests (Kaufman et al., 2013). This class extends Hoeffding's test by considering data-driven (or all) partitions of an arbitrary order \code{KxK} of the paired sample (or rank) sapce. Being distribution-free, it does not require resampling for computing p-values, making it considerably faster to apply than the HHG test (see example below).

Several different variants are implemented, which can be chosen by specifying the \code{variant} argument. The options are: 'spr.obs', 'spr.all', 'ppr.22.obs', 'ppr.33.obs', 'ppr.22.all', 'ppr.33.all', 'tpr.obs', 'tpr.all', 'ddp.obs', or 'ddp.all'. \emph{obs} variants implement DDP (data driven partition) tests, while \emph{all} variants implement ADP (all data partitions) tests. 

\emph{spr}, \emph{ppr}, and \emph{tpr} are optimized versions of the test for \code{K=2,3,4}, which are of \code{O(n^1), O(n^2), O(n^3)} time complexity, respectively (with \code{n = length(x)}). In the case of 'ppr' (\code{K=3}), there are two subvariants: '33' considers entire 3x3 contingency tables (similar to the test implemented for other \code{K} values); the '22' variant only looks at nonuniform distribution in the partition formed by the center cell versus the outer cells in every original 3x3 partition. The test for general \code{K} is implemented by the 'ddp' variants (so 'ddp.obs' is the general DDP test, and 'ddp.all' is the general ADP test). In this case \code{K} itself is given as a separate argument.
}

\value{
The optimized variants 'spr', 'ppr', and 'tpr' return four statistics:

\code{sum.chisq} - normalized sum of Pearson chi-squared statistics from the KxK contingency tables considered.

\code{sum.lr} - normalized sum of liklihood ratio ("G statistic") values from the KxK tables.

\code{max.chisq} - maximum Pearson chi-squared statistic from any of the KxK tables.

\code{max.lr} - maximum G statistic from any of the KxK tables.

Sum statistics are normalized by the total number of partitions multiplied by \code{length(x)}, which is the scale for measuring mutual information (in natural digits, aka 'nats').

As described in Kaufman et al. (2013), the \code{K > 4} variants 'ddp' and 'adp' resort to by-cell rather than by-partition enumeration, and thus do not currently return 'max' statistics (only 'sum' statistics are returned).
}

\references{
A Family of Consistent Univariate Distribution-Free Tests of Independence based on Sample Space Partitions. \emph{In preparation}.

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society Series} B, \bold{57}, 289-300.
}

\author{
  Shachar Kaufman.
}

\note{
The tests are designed for continuous random variables, so there should be no or very few ties in the input. Any ties in \code{x} and \code{y} are automatically resolved by random assignment of ranks (see the built-in function \code{rank}).
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
    ddp.nullsim[i] = xdp.test(x, y, variant = 'ppr.33.obs')$sum.lr
  }
  
  ddp.critical.value = quantile(ddp.nullsim, probs = 1 - bonferroni.level)
} else {
  # In the interest of saving time, this is what we get
  ddp.critical.value = 0.01569162
}

# Now we compute the observed test statistics for all pairs
compute.xdp.for.pair = function(pair.idxs) {
  x = as.numeric(hughes[pair.idxs[1], ])
  y = as.numeric(hughes[pair.idxs[2], ])
  return (xdp.test(x, y, variant = 'ppr.33.obs')$sum.lr)
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
