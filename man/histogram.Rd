\name{histogram}
\alias{histogram}
\encoding{latin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{histogram with automatic choice of bins}
\description{Construction of regular and irregular histograms with different options for 
choosing the number and widths of the bins. By default, both a regular and an irregular 
histogram using a data-dependent penalty as described in detail in Rozenholc/Mildenberger/Gather (2009) are constructed. 
The final estimate is the one with the larger penalized likelihood.
}
\usage{
histogram(y, type = "combined", grid = "data", 
breaks = NULL, penalty = "default", 
greedy = TRUE, right=TRUE, freq=FALSE, control = list(), 
verbose = TRUE, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of values for which the histogram is desired.  }
  \item{type}{use \code{"irregular"} for an irregular and \code{"regular"} 
  for a regular histogram. If \code{type="combined"} (default value) both a regular and an irregular histogram are computed 
  and the one with the larger penalized likelihood is chosen, see details below.} 
  \item{grid}{if \code{type="irregular"}, \code{grid} chooses the set of possible 
  partitions of the data range.
  The default value \code{"data"} gives a set of partitions constructed from the 
  data points, \code{"regular"} uses a fine regular grid of points as possible break points.  
  A regular quantile grid can be chosen using \code{"quantiles"}. Has no effect for regular histograms.}
  \item{breaks}{controls the maximum number of bins allowed in a regular histogram, or the size of 
  the finest grid in an irregular histogram when \code{grid} is set to \code{"regular"} or \code{"quantiles"}. 
  Usually not needed since the maximum bin number and the size of the finest grid are calculated by a formula depending 
  on the sample size \code{n}; the defaults for this can be changed using the parameters \code{g1}, 
  \code{g2} and \code{g3} in the \code{control} argument. Also see \code{maxbin} in the control argument which gives an absolute upper bound 
  bound on the number of bins if \code{type="regular"}.}
  \item{penalty}{controls which penalty is used. See description of penalties below.}
  \item{greedy}{logical; if \code{TRUE} and \code{type="irregular"}, a subgrid of the finest grid 
  is constructed by a greedy step to make the search procedure feasible. Has no effect for regular histograms.}
  \item{right}{logical; if \code{TRUE}, the histograms cells are right-closed (left open) intervals.}
  \item{freq}{logical; if \code{TRUE}, the y-axis gives counts in case of a regular histogram, otherwise density values are given. For irregular histograms, the y-axis always shows the density. Unlike \code{hist()}, defaults to \code{FALSE}.}
  \item{control}{list of additional control parameters. Meaning and default values depend on settings of \code{type}, 
  \code{penalty} and \code{grid}. See below.}
  \item{verbose}{logical; if \code{TRUE} (default), some information is given during histogram construction and the resulting \code{histogram} 
  object is printed.}
  \item{plot}{logical. If \code{TRUE} (default), the histogram is plotted.}
  \item{...}{further arguments and graphical parameters passed to \code{hist()}.}
}
\details{
The \code{histogram} procedure produces a histogram, i.e. a piecewise constant density estimate 
from a univariate real-valued sample stored in a vector \code{y}. Let \eqn{n} denote the length of \code{y}. 
The range of the data is partitioned into \eqn{D} intervals - called bins - and the density estimate on the \eqn{i}-th bin is 
given by \eqn{N_i/(n*w_i)} where \eqn{N_i} is the number of observations in the \eqn{i}-th bin and \eqn{w_i} is its width. The histogram thus defined is the maximum likelihood estimate among all densities that are piecewise constant w.r.t. this partition. 
The arguments of \code{histogram} given above determine the way the partition is chosen. In a regular histogram, the partition
consists of \eqn{D} bins of the same widths, and the histogram is determined by the choice of \eqn{D}. Strategies based on different
criteria can be chosen using the \code{penalty} option. The maximum number of bins can be controlled by either the \code{breaks} 
argument or the entries \code{g1}, \code{g2} and \code{g3} in the \code{control} argument. 

An irregular histogram allows for bins of different widths. In this case, not only the number \eqn{D} of bins but also the breakpoints
between the bins must be chosen. The set of allowed breakpoints is given by the finest partition selected using the \code{grid} argument. 
At the moment a finest regular grid is supported (\code{grid="regular"}) as well as grids with possible breakpoints either equal 
to the observations or between the observations (\code{grid="data"} and \code{between} in the \code{control} argument set to 
\code{FALSE} or \code{TRUE}, respectively). Setting \code{grid="quantiles"} gives a grid based on regular sample quantiles. 
If the \code{breaks} argument is \code{NULL}, \deqn{G(n)=g1*n^{g2}*(\log(n))^{g3}} controls the grid in the following way: the smallest 
allowed bin width in a \code{"data"} grid is \eqn{1/G(n)} times the sample range, while for \code{grid="regular"} and 
\code{grid="quantiles"} the finest grid has \code{floor(G(n))} bins. The parameters \code{g1}, \code{g2} and \code{g3} 
can be changed by modifying the corresponding components in the \code{control} argument. If \code{breaks} is a positive number, 
its integer part is used instead of \code{G(n)}. Different strategies for selection of \eqn{D}
and the bin boundaries can be chosen using the \code{penalty} option.

To reduce calculation time for irregular histograms, a subset of the breakpoints of the finest grid can be chosen by starting 
from a one-bin histogram and then subsequently finding the split of an existing bin that leads to the largest increase in the 
loglikelihood. The full optimization is then performed only over all partitions with endpoints from the subset
thus constructed. This is achieved by setting \code{greedy=TRUE}.  To reduce calculation time for regular histograms, the \code{maxbin}
parameter in the \code{control} argument gives an upper bound for the number of bins. The default value is 1000. 

Using \code{type="combined"} (the default value), both a regular and an irregular histogram are constructed using a penalized likelihood approach
and the one with the larger penalized likelihood is chosen. In this case, the regular histogram is always constructed using the \code{br} penalty. 
The \code{penalty} parameter and all other options control the construction of the irregular histogram. \code{penalty} 
must be equal to \code{"penA"}, \code{"penB"} or \code{"penR"}, since otherwise 
comparison of penalized likelihood values would not be meaningful.
}
\value{
an object of class "histogram" which is a list with the same components as in the \code{\link{hist}} command.
  }
\section{Penalties}{ 
Most settings of \code{penalty} lead to a penalized maximum likelihood histogram. For a sample of size \eqn{n} 
and a partition \eqn{J} that divides the sample range into \eqn{D} bins, define \eqn{N_i} as the number of observations in the \eqn{i}-th bin,
\eqn{i=1,...,D} and \eqn{w_i} as the width of the the \eqn{i}-th bin, \eqn{i=1,...,D}. In this section, 
the index in sums and products is always \eqn{i=1,\ldots,D}. For any partition \eqn{J}, 
and a fixed sample, the penalized loglikelihood is defined as 
\deqn{\sum N_i * \log(N_i/(n * w_i))-pen(J).}
The possible penalties are:
\describe{
\item{\code{penA}}{Penalty given in formula (5) in in Rozenholc, Mildenberger and Gather (2009):
\deqn{pen(J)=c \log {{n-1} \choose {D-1}} + \alpha(D-1) + ck\log(D) + 2\sqrt{c\alpha(D-1)(\log{{n-1} \choose {D-1}} +k \log D) } ,}
where the default values are \eqn{c=1}, \eqn{\alpha=0.5} and \eqn{k=2}. These can be changed using the \code{c}, 
\code{alpha} and \code{k} components of \code{control}.
}
\item{\code{penB}}{Simplified version of formula (5) in Rozenholc, Mildenberger and Gather (2009):
\deqn{pen(J)=c \log {{n-1} \choose {D-1}} + \alpha(D-1) + \log^{2.5} D,}
where the default values are \eqn{c=1} and \eqn{\alpha=1}. These can be changed using the \code{c} and
\code{alpha} components of \code{control}. Default penalty for irregular and combined histograms.
}
\item{\code{penR}}{Data-dependent penalty as given in formula (6) in Rozenholc, Mildenberger and Gather (2009):
\deqn{pen(J)=c \log {{n-1} \choose {D-1}} + (\alpha/n) \sum N_i/w_i + \log^{2.5} D,}
where the default values are \eqn{c=1} and \eqn{\alpha=0.5}. These can be changed using the \code{c} and
\code{alpha} components of \code{control}.  
}
\item{\code{aic}}{Akaike's Information Criterion (AIC). Defined by \eqn{pen(J)=\alpha*D}, where \eqn{\alpha} 
is 1 by default and may be changed using the \code{alpha} parameter in the \code{control} argument.}
\item{\code{bic}}{Bayesian Information Criterion (BIC). Defined by \eqn{pen(J)=\alpha*\log(n)*D}, where \eqn{\alpha} 
is 0.5 by default and may be changed using the \code{alpha} parameter in the \code{control} argument.}
\item{\code{nml}}{Normalized Maximum Likelihood. Formula is given in Davies, Gather, Nordman, Weinert (2009). 
Only available for regular histograms.}
\item{\code{br}}{Improved version of AIC for regular histograms as given in Birge and Rozenholc (2006). Defined as 
\eqn{pen(J)=D + \log^{2.5} (D)}. Default penalty for regular histograms, not available for irregular histograms.}
}
Some settings of \code{penalty} do not lead to maximization of a penalized likelihood but optimzation of different measures. These are:
\describe{
\item{\code{cv}}{Leave-p-out crossvalidation. Different variants can be chosen by setting the \code{cvformula} and 
\code{p} components in the \code{control} argument. \code{cvformula=1} and \code{cvformula=2} are available both for 
regular and irregular histograms. These are different versions of leave-p-out L2-crossvalidation, where choice of 
a partition is achieved by minimizing
\deqn{2\sum N_i/w_i - (n+1)\sum N_i^2/(n*w_i), }
or
\deqn{2(n-p) \sum N_i/w_i - (n-p+1) \sum N_i^2/w_i }
respectively, see formulas (11) and (12) in Celisse and Robin (2008). Since formula \code{1} does not depend on \eqn{p}, if 
the control parameter \code{p} is set to a value greater than 1 \code{cvformula} is set to \code{2}. 
Kullback-Leibler crossvalidation can be performed by setting \code{cvformula=3}. This is only available if \eqn{p=1} and 
\code{type="regular"}. The number of bins chosen is the maximizer of 
\deqn{\sum N_i \log(N_i-1) + n \log(D),}
see remark 2.3 in Hall and Hannan (1988).  
}
\item{\code{sc}}{Stochastic Complexity criterion, only available for regular histograms. Number of bins is chosen by maximizing
\deqn{\prod N_i! D^n(D-1)!/(D+n-1)!,} see formula (2.3) in Hall and Hannan (1988). 
}
\item{\code{mdl}}{Minimum Description Length criterion, only available for regular histograms. Number of bins is chosen by maximizing
\deqn{\sum (N_i-0.5)\log(N_i-0.5)-(n-0.5D)\log(n-0.5D)+n\log D -0.5D \log n ,} see formula (2.5) in Hall and Hannan (1988). }
}
}

\section{Control}{
The control parameter is a list with different components that affect the construction of the histogram. 
Meaning and default values depend on setting of the other parameters.
\describe{
\item{\code{alpha}}{Coefficient of the number of bins in penalties \code{penA}, \code{penB}, \code{aic}, \code{bic}. 
Coefficient of the data-driven part in the \code{penR} penalty. }
\item{\code{between}}{logical; if \code{TRUE} and \code{grid="data"}, possible bin ends are put between the observations, if
\code{FALSE} (default) they are placed at the observations}
\item{\code{c}}{Controls the weight of the penalty component that corrects for the multiplicity 
of partitions with the same number of bins in irregular histograms; only used in penalties \code{penA}, \code{penB} and
\code{penR}.}
\item{\code{cvformula}}{determines the type of crossvalidation to be performed. Can take the values 1, 2 and 3. 1 and 2 correspond 
to different versions of L2 crossvalidation, while \code{cvformula=3} performs Kullback-Leibler crossvalidation, which is at the moment 
only available for regular histograms. Note that \code{cvformula=3} automatically forces every bin to include at least 2 observations. 
If \code{p} is set to a value greater than 1, \code{cvformula=2} is used automatically.}
\item{\code{g1}}{The parameters \code{g1}, \code{g2} and \code{g3} control the maximum number of bins in a regular histogram as well as the
bin width and/or number for irregular histograms. Define \deqn{G(n)=g1*n^{g2}*(\log(n))^{g3}.} The maximum number  
of bins allowed in a regular histogram is given by \code{floor(G(n))}, the finest grid in an irregular histogram with \code{grid="regular"} is
obtained by dividing the sample range into \code{floor(G(n))} equisized bins, and if \code{grid="quantiles"}, the finest grid is 
obtained by dividing the interval \eqn{[0,1]} into equisized intervals and using the sample quantiles corresponding to the boundary points.
For an irregular histogram with \code{grid="data"}, a mimimum allowed bin size of \eqn{1/G(n)} is enforced. This can be disabled by
setting \code{g3} to \code{Inf}, causing \eqn{1/G(n)} to be zero. Default settings are \code{g1=1} and \code{g2=1} for all grids. 
Default values for \code{g3} are \code{-1} for \code{grid="regular"} and \code{grid="quantiles"} and \code{Inf} 
for \code{grid="data"}. Also see \code{maxbin}.}
\item{\code{g2}}{see \code{g1}. }
\item{\code{g3}}{see \code{g1}. }
\item{\code{k}}{Tuning parameter that only has an effect if \code{penalty="penA"}. Default value is 2. }
\item{\code{maxbin}}{Gives an absolute upper bound on the number of bins in order to keep the calculations feasible for large data sets. 
If the number of bins chosen via \code{breaks} or \code{g1}, \code{g2} and \code{g3} exceeds \code{maxbin}, \code{maxbin} in used as 
the maximum number of bins. Only has an effect for regular histograms. Defaults to 1000.}
\item{\code{p}}{Controls the number p of data points left out in the crossvalidation. Can take integer values between \code{1} 
(default) and \code{n-1}. If a value greater than 1 is chosen, \code{cvformula} is automatically set to 2 since crossvalidation 
formula 1 does not depend on p and Kullback-Leibler crossvalidation is only supported for \code{p=1}.}
\item{\code{quanttype}}{Determines the way the quantiles are calculated if \code{grid="quantiles"}. Corresponds to the \code{type}
argument in \code{\link{quantile}}, whose default \code{7} is also the default here.}
%\item{\code{mincount}}{integer between 0 and \code{floor(n/2)} that controls the minimum number of observations per bin. 
%Defaults to \code{0} \code{penalty=="mdl"} or \code{penalty=="cv"} and \code{cvformula==3} in the regular case, where the minimum required
%number of counts per bin is 1 and 2, respectively.}
}
}

\references{
Birg?, L. and Rozenholc, Y. (2006). How many bins should be put in a regular histogram? ESAIM: Probability and Statistics,
10, 24-45. 

Celisse, A. and Robin, S. (2008). Nonparametric density estimation by exact leave-p-out
cross-validation. Computational Statistics and Data Analysis 52, 2350-2368.

Davies, P. L., Gather, U., Nordman, D. J., and Weinert, H. (2009): A comparison of
automatic histogram constructions. ESAIM: Probability and Statistics, 13, 181-196.

Hall, P. and Hannan, E. J. (1988). On stochastic complexity and nonparametric density estimation. Biometrika 75, 705-714.

Rozenholc, Y, Mildenberger, T. and Gather, U. (2009). Combining regular and irregular histograms by penalized likelihood. 
Discussion Paper 31/2009, SFB 823, TU Dortmund.
\url{https://eldorado.tu-dortmund.de/handle/2003/26529}

Rozenholc, Y., Mildenberger, T., Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. Computational Statistics and Data Analysis 54, 3313-3323.
}
\author{Thoralf Mildenberger, Yves Rozenholc, David Zasada.}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{hist}}}
\examples{

## draw a histogram from a standard normal sample

y<-rnorm(100)
histogram(y)

## draw a histogram from a standard exponential sample

y<-rexp(1500)
histogram(y)

## draw a histogram from a normal mixture 

n<-sum(sample(c(0,1),1000,replace=TRUE))
y<-c(rnorm(n,mean=5,sd=0.1),rnorm(1000-n))
histogram(y)

## the same using a regular histogram with Kullback-Leibler CV 

n<-sum(sample(c(0,1),1000,replace=TRUE))
y<-c(rnorm(n,mean=5,sd=0.1),rnorm(1000-n))
histogram(y,type="regular",penalty="cv",control=list(cvformula=3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{smooth}% __ONLY ONE__ keyword per line
