% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PS_prior.R
\name{PS_prior}
\alias{PS_prior}
\alias{PS_prior.default}
\alias{PS_prior.beta}
\alias{PS_prior.norm}
\title{Calculating the Propensity Score-Integrated Informative Priors}
\usage{
PS_prior(
  formula,
  data,
  outcome,
  study,
  treat,
  method,
  distance,
  ratio,
  ps.method,
  trim
)

PS_prior.default(
  formula,
  data,
  outcome,
  study,
  treat,
  method,
  distance,
  ratio,
  ps.method,
  trim
)

PS_prior.beta(
  formula,
  data,
  outcome,
  study,
  treat,
  method,
  distance,
  ratio,
  ps.method,
  trim
)

PS_prior.norm(
  formula,
  data,
  outcome,
  study,
  treat,
  method,
  distance,
  ratio,
  ps.method,
  trim
)
}
\arguments{
\item{formula}{A two-sided \code{\link{formula}} object containing the study
indicator and covariates to be used in creating the distance measure used
in the matching. This formula will be supplied to the functions that estimate
the distance measure. For example, the formula should be specified as
\code{G ~ X1 + X2 + \dots} where \code{G} represents the name of study indicator
 and \code{X1} and \code{X2} are covariates.}

\item{data}{A data frame containing the variables named in \code{formula}
and possible other arguments.}

\item{outcome}{The variable name of the outcome.}

\item{study}{The variable name of the study indicator.}

\item{treat}{The variable name of the treatment indicator.}

\item{method}{The matching method to be used. The allowed methods are
\code{"nearest"} for nearest neighbor matching (on
the propensity score by default), \code{"optimal"} [\link[MatchIt]{method_optimal}]
for optimal pair matching, \code{"full"} [\link[MatchIt]{method_full}] for optimal
full matching, \code{"genetic"} [\link[MatchIt]{method_genetic}] for genetic
matching, \code{"cem"} [\link[MatchIt]{method_cem}] for coarsened exact matching,
\code{"exact"} [\link[MatchIt]{method_exact}] for exact matching,
\code{"cardinality"} [\link[MatchIt]{method_cardinality}] for cardinality and
template matching, and \code{"subclass"} [\link[MatchIt]{method_subclass}] for
subclassification. When set to \code{"NULL"}, no matching will occur, but
propensity score estimation and common support restrictions will still occur
if requested. See the linked pages for each method for more details on what
these methods do, how the arguments below are used by each on, and what
additional arguments are allowed.}

\item{distance}{The distance measure to be used. Can be either the name of a
method of estimating propensity scores (e.g., \code{"glm"}), the name of a
method of computing a distance matrix from the covariates (e.g.,
\code{"mahalanobis"}), a vector of already-computed distance measures, or a
matrix of pairwise distances. See [\link[MatchIt]{distance}] for allowable
options. The default is \code{"glm"} for propensity scores estimated with
logistic regression using \code{glm()}. Ignored for some methods; see individual
methods pages for information on whether and how the distance measure is
used.}

\item{ratio}{For methods that allow it, how many historical control units should be
matched to each current control unit in \eqn{k:1} matching. Should be a single integer
value. See the individual methods pages for information on whether and how
this argument is used. The default is 1 for 1:1 matching.}

\item{ps.method}{PS method utilize to calculate an informative prior
based on historical data. The allowed methods are \code{"Weighting"}
or \code{"Matching"}. The default method is \code{"Weighting"}.}

\item{trim}{Lower and upper bound of trimming used in \code{"Weighting"}.
The default is [0.1,0.9].}
}
\value{
Displays the informative prior calculated from historical data
based on the selected PS method.
}
\description{
The \code{PS_prior} function is designed to calculate the Propensity Score-Integrated
(PS) informative prior constructed based on historical data.
}
\details{
This function aims to calculate informative priors using historical
 data by incorporating covariate information to enhance borrowing strength
 and address prior-data conflicts.

 Let \eqn{G} be the study indicator, where \eqn{G = 1} indicate patient is
 from current control study, and \eqn{G = 0} indicate patient is from
 historical control study. Given the covariates data \eqn{X}, the propensity
 score is defined as follows,
 \deqn{e(X) = \Pr(G = 1 | X),}
 where \link[MatchIt]{distance} allows different methods to estimate the
 propensity scores.

 Calculate informative prior through PS matching is to identify a subset of
 historical data (\eqn{D_h^*}) that have similar PS as current control data
 (\eqn{D}). Various algorithms are available for PS matching, please refer
 to \code{method}. The informative prior can then be calculated based on
 the matched historical dataset.

 Alternative, we can utilize the inverse probability of treatment weighting
 (IPTW) to adjust the distribution of \eqn{X} in historical data \eqn{D_h},
 making it similar to that in \eqn{D}. Specifically, for the \eqn{i}th
 subject, we assign a weight \eqn{\alpha_i} to the outcome \eqn{y_i} in
 \eqn{D_h} based on its PS \eqn{e(X_i)} and a fixed weight \eqn{\alpha_i = 1}
 to \eqn{X_i} in \eqn{D}, as follows:
 \deqn{\alpha_i = G_1 + (1 - G_i) \frac{e(X_i)}{1 - e(X_i)}.}
 To avoid extremely large weights that may compromise IPTW, symmetric
 trimming rule can be used to trim the tails of the PS distribution by
 input \code{trim} with default [0.1,0.9], that is to trim observations
 whose estimated PS is outside of this range.

 To standardized \eqn{\alpha}, we compute the effective sample size (ESS),
 which approximately reflects the level of precision or equivalently its
 sample size, retained in the sample after weight as
 \eqn{n^{*}_h = (\sum \alpha_i)^2 / \sum{\alpha_i^2}}. The standardized weight
 is given by
 \deqn{\alpha_i^{*} = G_i + (1 - G_i)\frac{G_i}{\sum{\alpha_i} / n_h^{*}}.}

 For binary endpoint \eqn{Y \sim Ber(\theta)}, the informative
 prior \eqn{\pi_1(\theta)} can be constructed as follows,
 \deqn{\pi_1(\theta) \propto L(\theta | D_h, \alpha^{*}) \pi_0(\theta)
  = Beta(a + \sum \alpha_i^{*}y_i, b + n_h^* - \sum \alpha_i^{*}y_i )\},}
 where \eqn{\pi_0(\theta)} is a non-informative prior, a natural choice is
 \eqn{Beta(a, b)}, with \eqn{a = b = 1}.

 For continuous endpoint \eqn{Y \sim N(0, \sigma^2)}, suppose \eqn{\sigma^2}
 is unknown, with non-informative prior \eqn{p(\theta, \sigma^2) \propto 1/\sigma^2},
 \eqn{\pi_1(\theta)} follows a student-\eqn{t} distribution with degree of
 freedom \eqn{n_h^{*} - 1}. Given that \eqn{n_h^{*}} is moderate and large,
 it can be approximated by a normal distribution
 \eqn{N(\bar{y}^{*}, {s^{*}}^2 / n_h^{*})} with
 \deqn{\bar{y}^{*} = \sum \alpha_i^* y_i / \alpha_i^*, ~~ {s^{*}}^2 =
 \sum \alpha_i^* (y_i - \bar{y}^{*})^2 / (n_h^{*}  - 1).}
}
\section{Functions}{
\itemize{
\item \code{PS_prior.default()}: The function calculates the Propensity Score-Integrated
informative prior based on historical data for binary and continuous endpoint.

\item \code{PS_prior.beta()}: The function calculates the Propensity Score-Integrated
informative prior based on historical data for binary endpoint.

\item \code{PS_prior.norm()}: The function calculates the Propensity Score-Integrated
informative prior based on historical data for continuous endpoint.

}}
\examples{
## Load example data
data('PS_SAM_data')
## Subset the data to contain historical data and current control
dat <- PS_SAM_data[PS_SAM_data$A == 0, ]
str(dat)

## Examples for binary endpoints
## Generate the informative prior based on historical data using PS Matching
summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
                 data = dat, ps.method = 'Matching', method = 'nearest',
                 outcome = 'Y_binary', study = 'G', treat = 'A'))

## Generate the informative prior based on historical data using PS Weighting
summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
                 data = dat, ps.method = 'Weighting',
                 outcome = 'Y_binary', study = 'G', treat = 'A'))

## Examples for continuous endpoints
## Generate the informative prior based on historical data using PS Matching
summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
                 data = dat, ps.method = 'Matching', method = 'nearest',
                 outcome = 'Y_continuous', study = 'G', treat = 'A'))

## Generate the informative prior based on historical data using PS Weighting
summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
                 data = dat, ps.method = 'Weighting',
                 outcome = 'Y_continuous', study = 'G', treat = 'A'))


}
\references{
Zhao Y, Laird G, Chen J, Yuan Y. PS-SAM: doubly robust
propensity-score-integrated self-adapting mixture prior to dynamically
borrow information from historical data.
}
\seealso{
\link[MatchIt]{matchit}
}
