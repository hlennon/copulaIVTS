\name{kalman_filter}
\alias{kalman_filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kalman Filter}

\description{
Kalman Filter for fitting an arma(p,q) model to data
Estimated parameters perfectly well for large sample sizes
arima() function much quicker but doesn't give tau_t or e_t's though just final parameter estimates}
\usage{
kalman_filter(F, G, Q, R, Y, X0, P0)
}
\arguments{
  \item{F}{a numeric value }
  \item{G}{a numeric value }
  \item{Q}{a numeric value}
  \item{R}{a numeric value }
  \item{Y}{a numeric value}
  \item{X0}{Initial Observation Matrix}
  \item{P0}{a numeric value}
}
\details{
The Kalman_filter, to be used to calculate the arma likelihood}

\value{
Outputs a vector of e_t's and vector of tau_t's, both length n
}
\references{
https://www.ma.man.ac.uk/~jy/tsaf/
}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
This is the R version of Dr Jingsong Yuan's Code written in Matlab
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
kalman.m
}
\examples{
library(forecast)
data(ins.data)
x <- ins.data[1:500]
n <- length(x)

p <- 1
q <- 1

a <- 0.7
b <- 0.5

# Set up Kalman Filter to compute e_t's and tau_t's (aka S)
if (p < q + 1){                                              # Make sure p = q + 1
            a <- c(a, array(0,c(q + 1 - p)))                 # pad a with zeros
            p <- q + 1
            }                                                # increase p to q+1
if (p > q + 1){
            b <- c(b, array(0,c(p - 1 - q))) 	             # pad b with zeros
            q = p - 1;                	    	 # increase q to p-1
            }

r      <- ARMAacf(ar = a, ma = b, pacf = FALSE)
sigma2 <- 1/r[1]

F  <- rbind(a, cbind(diag(1,(p-1)), c(rep(0,(p-1)))))
G  <- rbind(c(1, b))
Q  <- diag(c(1, array(0,c(p-1, 1))))
R  <- 0
X0 <- array(0, c(p,1))
P0 <- toeplitz(r[1:p])

kal<- kalman_filter(F, G, Q, R, x, X0, P0)
attributes(kal)
}


