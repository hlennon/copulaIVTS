\name{MCEM}
\alias{MCEM}
\title{Monte Carlo EM algorithm
}
\description{
Use: Computes the MLEs of the Gaussian copula model using the Monte Carlo EM algorithm
}


\usage{
MCEM(obs.data, Monte_Carlo_samples, p, q, initial,  num_iterations,
                prec=0.01,  marginal_dist="negbin",
                optim_method="nmkb", compute_stderrors=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs.data}{ Observed Data }
  \item{Monte_Carlo_samples}{The number of Monte Carlo samples, m, to approximate the conditional expectation, if m is a vector it repeats num_iterations with each numbr in the vector}
  \item{p}{order of the ar parameters in ARMA(p,q)}
  \item{q}{order of the ma parameters in ARMA(p,q)}
  \item{initial}{ initial parameter values for the EM algorithm }
  \item{prec}{ precision for stopping criteria }
  \item{num_iterations}{ Number of iterations}
  \item{marginal_dist}{ Marginal distributions}
  \item{optim_method}{ Default is "nmkb" from the R package dfoptim, however "CG", "BFGS" or "nlm" can be used. optimx() from the optimx R package suggests all 4 can be used but in practice only "L-BFGS-B" have provided a good optimisation routine}
    \item{compute_stderrors}{TRUE or FALSE}
}
\details{
Default is "L-BFGS-B", however "CG", "BFGS" or "nlm" can be used. optimx() from the optimx R package suggests all 4 can be used but in practice only "L-BFGS-B" have provided a good optimisation routine.

There are two types of optimisation: function minimisation and mathematical programming. Most optimisers in R are function minimisers. method="L-BFGS-B" is the method of Bryd et al. (1995) which is a modification of BFGS quasi-newton method.
"Dela_lx" for change in log-likelihood of observed data, as following Chen and Ledolter (1995)

}
\value{
Iterations for the MCEM algorithm which converge to the MLE estimates of the ARMA(1,1) parameters
}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}

\references{
Lennon H., & Yuan J., Estimation of a digitised Gaussian ARMA model by Monte Carlo Expectation Maximisation, Computational Statistics & Data Analysis 2018;8:e020683

Lennon, H., 2016. Gaussian copula modelling for integer-valued time series (Doctoral Thesis, University of Manchester).
}
\note{
Produces pdf plots
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
First and Second Year Transfer Reports}
\examples{
library(numDeriv)
library(MASS)
library(VGAM)
library(dfoptim)
library(polynom)
library(dfoptim)
ls()

set.seed(1)
obs.data <- generate_data(10, 0.2, 0.1, "nbinom", c(3, 0.4))$x


}

