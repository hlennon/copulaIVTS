\name{arma_loglikelihood}
\alias{arma_loglikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{arma_loglikelihood
}
\description{
Evaluates the ARMA(p,q) log-likelihood function at the set of arma parameters
e.g.
order = c(1,1)
data = arma data
parameters = arma parameters (phi_1, ..., phi_p, theta_1, ..., theta_q)
}
\usage{
arma_loglikelihood(parameters, data, ord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{a vector of arma parameters}
  \item{data}{a vector of arma data}ord
 \item{ord}{ order of arma model}
}

\details{
Evaluates the ARMA(p,q) logliklihood function.
'data' and 'order' must be specified as global variables }

\value{
scalar value of the arma LOG-likelihood evaluated at vector (parameters) for the given data}


\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
Evaluates Reduced Gaussian log-likelihood
g   <-  n*log(2*pi) + n*log(sigma2) + sum(log(kal$S)) + sum(kal$E^2/kal$S)/sigma2
g   <- -g/2
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{arma_lik}

\examples{
library(forecast)

# Example 1
data(ins.data)
data  <- ins.data[1:200]
order <- c(1,1)
arma_loglikelihood(c(0.7, -0.5), data, c(1,1))



# Example 2
data  <- arima.sim(model=list(ar =c(0.4, -0.5), ma=c(0.8, 0.2)), n=2000)
order <- c(2,2)
para  <- c(0,0,0,0)

# takes 1 minute to optimise
# optim(para, arma_loglikelihood,  data=data, ord=order, control=list(fnscale=-1))$par
# [1]  0.4569347 -0.5244361  0.7507188  0.1448333



}
