\name{GHK}

\alias{GHK}

\title{GHK simulator}

\description{
 Simulates m samples of the n-dimensional zero mean Truncated Normal
 MVTN(mu=0, Sigma, a, b) via the Geweke-Hajivassiliou-Keane (GHK) sampler.
 The GHK simulation method is EXACT, in comparison to the rtmvnorm() function in the {tmvtnorm} R package which has the arguement algorithm taking values c("rejection", "gibbs", "gibbsR").
 The methods "gibbs" and "gibbsR" are very fast but only provide an approximate method, rather than an exact draw from an n-dimensional zero mean Truncated Normal
 MVTN(mu=0, Sigma, a, b).
###################################################################################
## GHK Simulator to exactly generate MVTN Samples
## Simulates n samples of the m-dimensional Truncated Normal
## MVTN( 0, Sigma, a, b)
## Zero Mean of MV Normal distn
## Sigma is the Covariance Matrix (m x m) of multivariate Normal distribution
## a is the lower truncation (m x 1 vector)
## b is the upper truncation (m x 1 vector)
##
## Arguements to be used:
## n - desired number of samples
## Sigma (Cov matrix)
## a (lower truncation) - which can be set to -Inf vector for right Truncation
## b (upper truncation) - which can be set to +Inf vector for left Truncation
##
## Requires: (no packages) m, Sigma, a, b,
## Output:   m random samples from n-dimensional TN(0, Sigma, a b) distribution.
}
\usage{
GHK(m, Sigma, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{      m is the number of samples generated }
  \item{mu=0}{   mu is set to be a vector of zeros of length n }
  \item{Sigma}{  Sigma is the Covariance Matrix (n x n) of multivariate Normal distribution }
  \item{a}{      a is the lower truncation (n x 1 vector)}
  \item{b}{      b is the upper truncation (n x 1 vector)}
}

\details{
The generation of random numbers from a truncated multivariate normal distribution is done using the GHK simulator. The m samples are output into m columns. Each column is one sample from the MV TN(mu=0, Sigma, a, b).


}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
http://www.hss.caltech.edu/~mshum/gradio/ghk_desc.pdf}
\author{
Hannah Lennon
}
\note{
The function rtmvnorm() is used from the mvtnorm R Package.}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# A 2-dimensional example
Sigma  <- toeplitz(c(1, 0.5))
a <- c(-1, -2)
b <- c( 1,  1)
GHK(3, Sigma, a, b)

# Test to
cbind(a, GHK(5, Sigma, c(-1, -2), c(1, 1)), b)

# compare to approximate draw
library(tmvtnorm)
t(rtmvnorm(3, mean=c(0,0), Sigma, a, b))


# A 300-dimensional example
d = 300
Sigma  <- toeplitz(c(0.9^(1:d)))
a <- rep(-2, times=d)
b <- rep( 1, times=d)
system.time(GHK(d, Sigma, a, b))

# compare to approximate draw
library(tmvtnorm)
system.time(t(rtmvnorm(d, mean=rep(0,times=d), Sigma, a, b ,algorithm="gibbsR")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
