% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAM_weight.R
\name{SAM_weight}
\alias{SAM_weight}
\alias{SAM_weight.betaMix}
\alias{SAM_weight.normMix}
\alias{SAM_weight.gammaMix}
\title{Calculating Mixture Weight of SAM Priors}
\usage{
SAM_weight(if.prior, theta.h, method.w, prior.odds, data, delta, ...)

\method{SAM_weight}{betaMix}(if.prior, theta.h, method.w, prior.odds, data, delta, n, r, ...)

\method{SAM_weight}{normMix}(
  if.prior,
  theta.h,
  method.w,
  prior.odds,
  data,
  delta,
  m,
  n,
  sigma,
  ...
)

\method{SAM_weight}{gammaMix}(if.prior, theta.h, method.w, prior.odds, data, delta, u, w, ...)
}
\arguments{
\item{if.prior}{Informative prior constructed based on historical data,
represented (approximately) as a mixture of conjugate distributions.}

\item{theta.h}{Estimate of the treatment effect based on historical data.
If missing, the default value is set to be the posterior mean estimate from
\code{if.prior}.}

\item{method.w}{Methods used to determine the mixture weight for SAM priors.
The default method is "LRT" (Likelihood Ratio Test), the alternative option
is "PPR" (Posterior Probability Ratio). See Details section for more information.}

\item{prior.odds}{The prior probability of \eqn{H_0} being true compared to
the prior probability of \eqn{H_1} being true using PPR method. The default
value is 1. See Details section for more information.}

\item{data}{Data of the control arm from the current trial, see Methods
section for more details.}

\item{delta}{Clinically significant difference used for the SAM prior.}

\item{...}{Additional parameters required for different endpoints.}

\item{n}{Number of subjects in the control arm for continuous endpoint.}

\item{r}{Number of responses in the control arm for binary endpoint.}

\item{m}{Mean estimate in the control arm for continuous endpoint.}

\item{sigma}{Standard deviation in the control arm for continuous endpoint.}

\item{u}{Number of events in the control arm for time-to-event endpoint.}

\item{w}{Total observed time in the control arm for time-to-event endpoint.}
}
\value{
The mixture weight of the SAM priors.
}
\description{
The \code{SAM_weight} function is designed to calculate the mixture
weight of the SAM priors according to the degree of prior-data
conflicts (\emph{Yang, et al., 2023}).
}
\details{
SAM prior is constructed by mixing an informative prior
\eqn{\pi_1(\theta)}, constructed based on historical data, with a
non-informative prior \eqn{\pi_0(\theta)} using the mixture weight
\eqn{w} determined by \code{\link{SAM_weight}} function to achieve the
degree of prior-data conflict (\emph{Schmidli et al., 2015, Yang et al., 2023}).

 Let \eqn{\theta} and \eqn{\theta_h} denote the treatment effects
 associated with the current arm data \eqn{D} and historical data \eqn{D_h},
 respectively. Let \eqn{\delta} denote the clinically significant difference
 such that if \eqn{|\theta_h - \theta| \ge \delta}, then \eqn{\theta_h} is
 regarded as clinically distinct from \eqn{\theta}, and it is therefore
 inappropriate to borrow any information from \eqn{D_h}. Consider two
 hypotheses:

 \deqn{H_0: \theta = \theta_h, ~ H_1: \theta = \theta_h + \delta ~ or ~  \theta = \theta_h - \delta.}
 \eqn{H_0} represents that \eqn{D_h} and \eqn{D} are consistent (i.e.,
 no prior-data conflict) and thus information borrowing is desirable,
 whereas \eqn{H_1} represents that the treatment effect of \eqn{D}
 differs from \eqn{D_h} to such a degree that no information should be
 borrowed.

 The SAM prior uses the likelihood ratio test (LRT) statistics \eqn{R} to
 quantify the degree of prior-data conflict and determine the extent of
 information borrowing.

 \deqn{R = P(D | H_0, \theta_h) / P(D | H_1, \theta_h) = P(D | \theta = \theta_h) / \max(P(D | \theta = \theta_h + \delta), P(D | \theta = \theta_h - \delta)) ,}
 where \eqn{P(D | \cdot)} denotes the likelihood function. An alternative
 Bayesian choice is the posterior probability ratio (PPR):

 \deqn{R = P(D | H_0, \theta_h) / P(D | H_1, \theta_h) = P(H_0) / P( H_1) \times BF, }
 where \eqn{P(H_0)} and \eqn{P(H_1)} is the prior probabilities of \eqn{H_0}
 and \eqn{H_1} being true. \eqn{BF} is the Bayes Factor that in this case
 is the same as the LRT.

 The SAM prior, denoted as \eqn{\pi_{sam}(\theta)}, is then defined
 as a mixture of an informative prior \eqn{\pi_1(\theta)}, constructed
 based on \eqn{D_h} and a non-informative prior \eqn{\pi_0(\theta)}:

 \deqn{\pi_{sam}(\theta) = w\pi_1(\theta) + (1-w)\pi_0(\theta),}
 where the mixture weight \eqn{w} is calculated as:

 \deqn{w = R / (1 + R).}

 As the level of prior-data conflict increases, the likelihood ratio
 \eqn{R} decreases, resulting in a decrease in the weight \eqn{w}
 assigned to the informative prior and thus a decrease in information
 borrowing. As a result, \eqn{\pi_{sam}(\theta)} is data-driven and
 has the ability to self-adapt the information borrowing based on the
 degree of prior-data conflict.
}
\section{Methods (by class)}{
\itemize{
\item \code{SAM_weight(betaMix)}: The function calculates the mixture weight of SAM
priors for beta mixture distribution. The input \code{data} can be patient-level
data (i.e., a vector of 0 and 1 representing the response status of each
patient) or summary statistics (i.e., the number of patients and the number
of responses).

\item \code{SAM_weight(normMix)}: The function calculates the mixture weight of SAM
priors for normal mixture distribution. The input \code{data} should be
a vector of patient-level observations. The input \code{data} can be
patient-level data (i.e., a vector of continuous response of each
patient) or summary statistics (i.e., the mean estimate, number of subjects,
and the standard deviation in the control arm).

\item \code{SAM_weight(gammaMix)}: The function calculates the mixture weight of SAM
priors for gamma mixture distribution. The input \code{data} can be
patient-level data (i.e., a matrix with the first row as the censoring
indicator and the second row recording the observed time) or summary
statistics (i.e., the number of uncensored observations \code{u} and
total observed time \code{w}).

}}
\examples{
set.seed(123)
## Examples for binary endpoints
## Example 1: no prior-data conflict
## Suppose that the informative prior constructed based on historical data is
## beta(40, 60)
prior.historical <- mixbeta(c(1, 40, 60))
## Data of control arm
data.control     <- rbinom(60, size = 1, prob = 0.42)
## Calculate the mixture weight of the SAM prior
wSAM <- SAM_weight(if.prior = prior.historical,
                   delta = 0.15,        ## Clinically significant difference
                   data = data.control  ## Control arm data
                   )
print(wSAM)

## Example 2: in the presence of prior-data conflict, where the current data
## has 12 responses in 60 patients
wSAM <- SAM_weight(if.prior = prior.historical,
                   delta = 0.15,    ## Clinically significant difference
                   ## Methods to determine mixture weight for the SAM priors
                   ## by Posterior Probability Ratio
                   method.w = 'PPR',
                   ## Prior odds of favoring no prior-data conflicts to
                   ## the presence of prior-data conflict
                   prior.odd = 1/9,
                   n = 60,          ## Number of patients in the control arm
                   r = 12           ## Number of responses in the control arm
                   )
print(wSAM)

## Example 3: in the presence of prior-data conflict, where the current data
## has 12 responses in 60 patients
wSAM <- SAM_weight(if.prior = prior.historical,
                   delta = 0.15, ## Clinically significant difference
                   n = 60,       ## Number of patients in the control arm
                   r = 12        ## Number of responses in the control arm
                   )
print(wSAM)

## Examples for continuous endpoints
## Example 1: no prior-data conflict
## Suppose that the informative prior constructed from historical data is
## N(0, 3)
sigma      <- 3
prior.mean <- 0
prior.se   <- sigma/sqrt(100)
prior.historical <- mixnorm(c(1, prior.mean, prior.se), sigma = sigma)
## Data of the control arm
data.control     <- rnorm(80, mean = 0, sd = sigma)
wSAM <- SAM_weight(if.prior = prior.historical,
                   delta = 0.3 * sigma,    ## Clinically significant difference
                   data = data.control     ## Control arm data
                   )
print(wSAM)

## Example 2: in the presence of prior-data conflict, where the current data
## has mean of 0.5
data.control     <- rnorm(80, mean = 1, sd = sigma)
wSAM  <- SAM_weight(if.prior = prior.historical,
                    delta = 0.3 * sigma,    ## Clinically significant difference
                    data = data.control     ## Control arm data
                    )
print(wSAM)

## Examples for survival endpoints
## Example 1: no prior-data conflict
## Suppose the survival times from historical data follows exp(1) distribution
## with random censoring time follows U(0.5, 5) distribution
T_hi <- rexp(100, rate = 1)
C_hi <- runif(100, min = 0.5, max = 5)
## Indicators of the uncensored events
delta_hi <- as.numeric(T_hi < C_hi)
## Observed survival times from historical data
U_hi     <- T_hi
U_hi[delta_hi == 0] <- C_hi[delta_hi == 0]
## Construct the informative prior based on simulated historical data
prior.historical <- mixgamma(c(1, sum(delta_hi), sum(U_hi)),
                             param = 'ab', likelihood = 'exp')
## Suppose the survival times from control data follows exp(0.95) distribution
## with random censoring time follows U(0.5, 5) distribution
T_ci <- rexp(100, rate = 0.95)
C_ci <- runif(100, min = 0.5, max = 5)
## Indicators of the uncensored events
delta_ci <- as.numeric(T_ci < C_ci)
## Observed survival times from control data
U_ci     <- T_ci
U_ci[delta_ci == 0] <- C_ci[delta_ci == 0]
## Data of the control arm
data.control     <- rbind(sum(delta_ci), sum(U_ci))
wSAM <- SAM_weight(if.prior = prior.historical,
                   delta = 0.2,            ## Clinically significant difference
                   data = data.control     ## Control arm data
                   )
print(wSAM)

## Example 2: in the presence of prior-data conflict, where the current survival
## times follows exp(2) distribution with random censoring time follows U(0.5, 5)
## distribution
T_ci <- rexp(100, rate = 2)
C_ci <- runif(100, min = 0.5, max = 5)
## Indicators of the uncensored events
delta_ci <- as.numeric(T_ci < C_ci)
## Observed survival times from control data
U_ci     <- T_ci
U_ci[delta_ci == 0] <- C_ci[delta_ci == 0]
## Data of the control arm
data.control     <- rbind(sum(delta_ci), sum(U_ci))
wSAM  <- SAM_weight(if.prior = prior.historical,
                    delta = 0.2,            ## Clinically significant difference
                    data = data.control     ## Control arm data
                    )
print(wSAM)

}
\references{
Yang P, Zhao Y, Nie L, Vallejo J, Yuan Y.
SAM: Self-adapting mixture prior to dynamically borrow information from
historical data in clinical trials. \emph{Biometrics} 2023; 79(4), 2857-2868.
}
