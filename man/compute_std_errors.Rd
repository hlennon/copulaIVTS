\name{compute_std_errors}
\alias{compute_std_errors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compute_std_errors
}
\description{
Computation of the standard errors for the MCEM algorithm for the GCNB model
}
\usage{
compute_std_errors(theta, obs.data,  m, p, q, marginals)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{a vector of MLE parameter estimates}
  \item{obs.data}{a vector of observed data}
  \item{m}{The number of monte carlo simulations required}
  \item{p}{The order of the ARMA(p,q) model}
  \item{q}{The order of the ARMA(p,q) model}
    \item{marginals}{how are the marginals estimated }

}

\details{
The standard errors of the MLE's are not automatically produced for the EM algorithm.
Equation (3) of Chen and Ledolter, (1995) is used to compute the errors.
The Hessian and the score of the ARMA(p,q) likelihood is required. The analytical solution is sought, however for the tie being a numerical approximation is used, via the R package (numDeriv).
The expectation is found by using a Monte Carlo average over m samples.
}

\value{
The standard errors of the MLE's of the ARMA(p,q) parameters for the GCNB model
}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
The standard errors for the true unobserved time series data (we know it because we simulated it), will be smaller than the standard error of the digitised time series
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{Lennon, H (2014)}
\references{Chen and Ledolter, (1995), Monte Carlo EM Estimation.}

\examples{
library(numDeriv)
library(forecast)


set.seed(1)
IV.data   <- generate_data(200, 0.7, -0.5, "negbin",  c(3, 0.3))
obs.data  <- IV.data$x

order <- c(1, 1)
theta <- arima(IV.data$y, order=c(1,0,1), include.mean=FALSE)$coef
m     <- 10
n     <- length(obs.data)

# Standard errors for the MLE's obtained via a MCEM algorithm
compute_std_errors(theta, obs.data,  m, 1, 1)


# Standard errors for the true unobserved data (we know it because we simulated it)
arima(IV.data$y, order=c(1,0,1), include.mean=FALSE)
}

