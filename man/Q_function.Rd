\name{Q_function}
\alias{Q_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Q(theta, theta*) objective function in the EM algorithm}

\description{ Evaluates the Q function in EM algorithm for the Gaussian copula model with NB marginals involving the terms E[YiYj | X=x] 

Computes the Conditional likelihood given the observed data
Computes e2/t by simulating yy (EXACTLY GHK) and averaging
Arguements to be used: para (true parameter) & para.star  (current parameter)

}


\usage{ Q_function(para, para.star, m, low, up)}

\arguments{
  \item{para}{ true parameters, theta, of the function}
  \item{para.star}{ current parameters, theta star, of the function. Note theta* are the current values of the parameters in the iterative algorithm and therefore are fixed values once specified}
  \item{low}{ lower bound of values, dependent only on data and marginal distribution the Negative Binomial}
  \item{up}{ upper bound of values, dependent only on data and marginal distribution the Negative Binomial}

 \item{m}{The number of Monte Carlo samples, m, to approximate the conditional expectation}
 }
\details{
Objective Q funciton of the MCEM Algorithm 

 p - order of the ar parameters in ARMA(p,q), FIXED TO p=1
 q - order of the ma parameters in ARMA(p,q), FIXED TO q=1
 
 }
\value{
A numeric value of the conditional expectation of the complete likelihood function given the observed data
}
\references{
Dempster et al, (1977)}
\author{
Hannah Lennon}
\note{
Requires:
LD,
ensure_stationarity_causality, 
GHK
}

\examples{
library(polynom)

para      <- c(0.9, 0.5)
para.star <- c(0.8, 0.5)

low  <- c(rep(-1, 10))
up   <- c(rep( 1, 10))

Q_function(para, para.star, m=1, low, up)




## The function is currently defined as

Q_function <-function(para, para.star, m, low, up){       
            p          <- 1
            q          <- 1
            N          <- length(low)
            para.star[1:p] <- ensure_causality_invertibility(para.star[1:p])         # Ensures parameters and current pa
            a.star     <- para.star[1:p]                      # ar current parameters
            b.star     <- para.star[(p+1):(p+q)]              # ma current parameters          
            r.star     <- as.numeric(ARMAacf(ar=a.star, ma=b.star, lag.max=N, pacf=FALSE))   #r.star <- arma_cov(a.star, b.star, N); r.star <- r.star/r.star[1]
            P          <- toeplitz(r.star[1:N])                # Covariance matrix for Mult Truncated Normal
            para[1:p]  <- ensure_causality_invertibility(para[1:p]) 
            a          <- para[1:p]
            b          <- para[(p+1):(p+q)]         
            r          <- as.numeric(ARMAacf(ar=a, ma=b, lag.max=N, pacf=FALSE))                     # Theoretical autocovariances for parameters
            output     <- LD(r,N-1)
            PHI        <- t(output[, 1:N])
            tau        <- output[,(N+1)]
            sim.data   <- GHK(m, P, low, up)
            xxx       <<- sim.data            
            YY         <- sim.data *t(sim.data)/m
            tau.mat    <- diag(1/tau)                           # Matix with tau^-1 in diagonal elements
            A          <- PHI *tau.mat*t(PHI)                # Matrix form of equations
            like       <- sum(A*YY) + sum(log(tau)) 
            return(like)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
