\name{direct_GCNB_ARMA_likelihood}
\alias{direct_GCNB_ARMA_likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{direct_GCNB_ARMA_likelihood
}
\description{Direct Likelihood evaluation of GCNB model with ARMA(p,q) correlation matrix
Requires: mnormt, mvtnorm
}

\usage{
direct_GCNB_ARMA_likelihood(covariances, nb.p, nb.s, data, method="mnormt")
}
\arguments{
  \item{covariances}{Theoretical ARMA(p,q) covariances corresponding to the dependence structure}
  \item{nb.p}{ X_t ~ NB (p, s)}
  \item{nb.s}{ X_t ~ NB (p, s)}
  \item{data}{ observed data}
  \item{method}{ a string of value "mnormt" or "mvtnorm" only}
  }
\details{
Evaluated likelihood function directly}
\value{
A numeric value of the evaluated likelihood function given the observed data and dependence parameters, and marginal parameters nb.p and nb.s.
}
\references{
Shi and Valdez, (2012), Longitudinal Modelling of Insurance Claims using Jittering.}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
ARMA(p,q) dependence structure is used.

Cannot be optimised using optim(numeric(n), fn=direct_GCNB_ARMA_likelihood, nb.p=0.3, nb.s=3, data=data) which i think is due to the sadvm() function it calls from Fortan.
direct_GCNB_AR1_likelihood can be optimised over one parameter but is SLOW for n > 8.}


\examples{

library(mnormt)
library(mvtnorm)


n <- 4
set.seed(n)
simdata    <- arima.sim(n, model=list(ar=0.7))
data       <- qnbinom(pnorm(simdata), prob=0.3, size=3)
cov        <- ARMAacf(ar=0.7, lag.max=n-1, pacf=FALSE)

direct_GCNB_ARMA_likelihood(cov, 0.3, 3, data, method="mnormt")
direct_GCNB_AR1_likelihood(0.7,  0.3, 3, data, method="mnormt")

direct_GCNB_ARMA_likelihood(cov,  0.3, 3, data, method="mvtnorm")
direct_GCNB_AR1_likelihood(0.7,  0.3, 3, data, method="mvtnorm")

optimise(f=direct_GCNB_AR1_likelihood, interval=c(-0.9, 0.9), nb.p=0.3, nb.s=3, data=data, method="mnormt",  maximum=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
