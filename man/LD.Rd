\name{LD}
\alias{LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Levinson Durbin algorithm
}
\description{
Levinson Dubin Algorithm, following Brockwell & Davis 1987, to compute phi_k's and arrange in phi.mat 
Also computes phi_kk's and tau's needed for likelihood.
}
\usage{
LD(r, n-1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a vector of autocovariances }
    \item{n}{a numeric to specify the length of vector taus/dimension of matrix phi_kks}
}
\details{
See Brockwell and Davis}
\value{
\item{phi.mat }{A n x n matrix of phi_kk's}
\item{s}{ a vector of tau's, of length n}
%% ...
}
\references{
Brockwell and Davis (1987)}
\author{
Hannah Lennon
}

\examples{
r      <- as.numeric(ARMAacf(ar=0.7, ma=-0.5, lag.max=10, pacf=FALSE))   
n      <- 5 
output <- LD(r, n-1)

PHI    <- t(output[, 1:n])    # matrix of n x n 
tau    <- output[,(n+1)]      # vector of n
PHI; tau
            
}
