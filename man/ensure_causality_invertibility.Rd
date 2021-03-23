\name{ensure_causality_invertibility}
\alias{ensure_causality_invertibility}
\title{ensure_causality_invertibility}



\description{
Test to ensure AR parameters are causal and MA parameters are invertible.
If the AR parameters are used, the function converts them to casual
If the MA parameters are used, the function converts them to invertible parameters
If ar parameters,  (phi_1, phi_2, ..., phi_p)
If ma parameter,  (theta_1, theta_2, ..., theta_q)
}

\usage{
ensure_causality_invertibility(parameters)
}

\arguments{
\item{parameters}{vector of AR or MA parameters, in the order of  }
}
\details{
Test to ensure AR parameters are causal and MA parameters are invertible
It converts them to casual/incertible  parameters if not
}

\value{
parameter values that are causal or invertible, ensureing unique solutions.
}

\references{
Jingsong Yuan, Time Series Lecture Notes,
https://www.ma.man.ac.uk/~jy/tsaf/
}


\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}

\seealso{
Brockwell and Davis (1987)
}



\examples{
library(polynom)

ensure_causality_invertibility(c(0.5, 1.3))

p    <- 3
q    <- 3
para <- c(0.2, 1, 0.2, 1, -0.5, 0.3)

para[1:p] <- ensure_causality_invertibility(para[1:p])
para[(p + 1):(p + q)] <- -ensure_causality_invertibility(-para[(p + 1):(p + q)])

## The function is currently defined as
function (parameters)
{
    L <- length(parameters)
    ps <- polynomial(c(-parameters[L:1], 1))
    ps <- solve(ps)
    if (all(abs(ps) < 1))
        return(parameters)
    else {
        j <- numeric(L)
        k <- 0
        vec <- numeric(L)
        ps[abs(ps) >= 1] <- 1/ps[abs(ps) >= 1]
        while (k < L) {
            j <- add_1(j)
            k <- sum(j)
            vec[k] <- vec[k] + (-1)^(k + 1) * prod(ps^j)
                      }
        return(as.numeric(vec))
           }
}


}

