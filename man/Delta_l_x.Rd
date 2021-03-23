\name{Delta_l_x}
\alias{Delta_l_x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delta_l_x, computes the approximate change in observed likelihood }
\description{
computes the approximate change in observed likelihood via the kalman filter
}
\usage{
Delta_l_x(xxx, para, para_old, order)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxx}{a matrix of simulated data, reused from the E-step to calc change in l_x }
    \item{para}{updated parameter }
      \item{para_old}{parameters from previous iteration }
        \item{order}{a vector of p, q }
   
}
\details{
To compute xxx, you should run the Q_function, 
xxx is a by-product assigned as a global variable because Q_function must return a scalar to be optimised via the optim() function


}
\value{
a numeric value of the change in the observed-data log likelihood
%% ...
}
\references{
Brockwell and Davis (1987)}
\author{
Hannah Lennon
}

\examples{
library(polynom)
para      <- c(0.9, -0.5)
para.star <- c(0.8, -0.5)
low  <- c(rep(-1, 10))
up   <- c(rep( 1, 10))
p    <- 1
q    <- 1
Q_function(para, para.star, m=10, low, up, p, q)  
## output is xxx 
# The evaluation of the Q_function (above) is to compute xxx, 
#     it is a by-product assigned as a global variable

Delta_l_x(xxx, c(0.9, -0.5), c(0.89, -0.48), c(1,1))
Delta_l_x(xxx, c(0.9, -0.5), c(0.9, -0.5), c(1,1))
}
