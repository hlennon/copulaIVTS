\name{arma_cov}
\alias{arma_cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical Autocovariances for ARMA(p,q) model, given by 
model: x(t)-a(1)x(t-1)-...-a(p)x(t-p)=e(t)+b(1)e(t-1)+...+b(q)e(t-q)
}
\description{
Computes the theoretical Autocovariances for parameters a,b of length n where 
a = ar parameters
b = ma parameters
Note: This is Dr Yuan's Matlab Code written in R 
}


\usage{
arma_cov(a, b, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ AR coefficients in a (0 if MA)}
  \item{b}{ MA coefficients in b (0 if AR)}
  \item{n}{n number of autocovariances required including variance }
}

\details{Requires: a,b,n, (could use p and q)

model: x(t)-a(1)x(t-1)-...-a(p)x(t-p)=e(t)+b(1)e(t-1)+...+b(q)e(t-q)
}

\value{
r is a vector of theoretical autocovariances of length n, 
output: variance and autocovariances of x(t) up to lag n-1 assuming variance of e(t) = 1

}

\author{
Hannah Lennon and Jingsong Yuan
}
\note{
Note: This is Dr Yuan's Matlab Code written in R 

The methods used follow Brockwell & Davis (1991, section 3.3). Their equations (3.3.8) are solved for the autocovariances at lags 0, …, max(p, q+1), and the remaining autocorrelations are given by a recursive filter.
 Note 1:
 arma.cov() assumes sigma^2 = 1.... it doesnot
 therefre arma.cov outputs sigma^2 (R(0), R(1), .. R(n-1))
 so we must divide by sigma^2 
 divide by R(0) to get the theoretical autocovariances R(0), R(1), etc
 r <- arma.cov(0.5, 0.1, 5)
 r <- r/r[1]

 Note 2: ARMAacf is much faster and gives the same result
 r <- arma.cov(0.7, -0.5, 10) ; r <- r/r[1]
 r <- as.numeric(ARMAacf(0.7, -0.5, lag=10))
 aram.cov (Yuans modified) is 2.5 times faster than arma.covmatlab 
 but hard to see the difference for n<1000
 
 }
\seealso{
ARMAacf gives a vector of autocorrelations, named by the lags. 
}
\examples{
arma_cov(0.7, -0.5, 5)
r <- arma_cov(0.7, -0.5, 5)
r <- r/r[1]


ARMAacf(ar=0.7, ma=-0.5, lag.max=5, pacf=FALSE)
## This gives the theoretical autoCORRELATION funcition (e.g. divide by Var(Y_t)=r[1])


## The function is currently defined as
function (a, b, n) 
{
    p <- length(a)
    q <- length(b)
    c <- ARMAtoMA(a, b, q)
    aa <- c(1, -a)
    A <- matrix(0, p + 1, p + 1)
    for (i in 0:p) {
        for (j in 0:p) {
            k <- abs(i - j)
            A[i + 1, k + 1] <- A[i + 1, k + 1] + aa[j + 1]
        }
    }
    b <- c(1, b)
    c <- c(1, c)
    d <- NULL
    for (i in 0:p) {
        d <- rbind(d, sum(b * c))
        b <- c(b[-1], 0)
    }
    r <- solve(A, d)
    for (j in (p + 1):n) {
        dum <- 0
        for (i in 1:p) {
            dum <- dum + a[i] * r[j - i + 1]
        }
        r <- c(r, dum + sum(b * c))
        b <- c(b[-1], 0)
    }
    return(r[1:(n + 1)])
  }
}