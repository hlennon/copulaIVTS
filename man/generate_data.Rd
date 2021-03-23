\name{generate_data}
\alias{generate_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{generate_data}
\description{generates data from the Gaussian copula model with Negative Binomial or Poisson marginals
of length n}
\usage{
generate_data(n, a, b, marginal_dist, marginal_parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{length of times series}
\item{a}{autoregressive parameters}
\item{b}{moving average parameters }
\item{marginal_dist}{"poisson" or "nbinom" or "zeroinf_poisson" or "zeroinf_nbinom"}
\item{marginal_parameters}{If negbin then c(nb.s, nb.p), If poisson then lambda}
}
\details{
Simulates 1 time series from the GCNB or GCP model
}
\value{
An integer-valued time series of length n
}
\references{
Schmidt (2007)
}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
Marginals fixed to be Poisson, nb, zero_inf_poisson or zero_inf_negbin at present.
The example below has ARMA(2,2) dependence structure, NB(p,s) marginal structure which is coupled together using a Gaussian copula
Note:

"nbinom" uses p and s representation.
"zeroinf_nbinom" uses negbinomial() which uses the alternate parametrization via mean mu and an index parameter k
An alternative parametrization (often used in ecology) is by the mean mu, and size, the dispersion parameter, where prob = size/(size+mu). The variance is mu + mu^2/size in this parametrization.




}


\examples{
library("VGAM")
generate_data(4, 0.7, -0.5, "poisson", 3)
generate_data(4, 0.7, -0.5, "nbinom",  c(3, 0.3))
generate_data(4, 0.7, -0.5, "zeroinf_poisson",  3)
generate_data(4, 0.7, -0.5, "zeroinf_nbinom",  c(3, 0.3))


tabulate(generate_data(40, 0.7, -0.5, "nbinom",  c(3, 0.5)) $x);
tabulate(generate_data(4, 0.7, -0.5, "zeroinf_nbinom",  c(3, 0.5)) $x)


set.seed(12)
data        <- generate_data(10000, c(0.2, 0.3), c(-0.1, 0.5), "nbinom",  c(5, 0.5))
ins.data21  <- data$x
ins.data21y <- data$y
hist(ins.data21)

## Code currentlt is defined as:
function (n, a, b, marginal_dist, marginal_parameters){

            if(marginal_dist=="poisson"){        lambda  <- marginal_parameters[1]}
            if(marginal_dist=="zeroinf_poisson"){lambda  <- marginal_parameters[1]}
            if(marginal_dist=="nbinom"){          nb.s   <- marginal_parameters[1]
                                              nb.p   <- marginal_parameters[2]}
            if(marginal_dist=="zeroinf_nbinom"){  nb.s   <- marginal_parameters[1]
                                                  nb.p   <- marginal_parameters[2]}


            ar.sigma <- arma_cov(a, b, 0)
            ar.sigma <- sqrt(1/ar.sigma)

            if(a[1]==0){  y <- arima.sim(list(ma = b), n, sd = ar.sigma)}
            else if(b[1]==0){  y <- arima.sim(list(ar = a), n, sd = ar.sigma)}
            else{ y <- arima.sim(list(ar = a, ma = b), n, sd = ar.sigma)}

            if(marginal_dist=="poisson")  x <- qpois(pnorm(y), lambda)
            if(marginal_dist=="nbinom")   x <- qnbinom(pnorm(y), size = nb.s, prob = nb.p)
            if(marginal_dist=="zeroinf_poisson") {
                        library(VGAM)
                        x <- qzipois(pnorm(y), lambda, pstr0 = 0)
            }
            if(marginal_dist=="zeroinf_nbinom") {
                        library(VGAM)
                        x <- qzinegbin(pnorm(y), size=nb.s, prob=nb.p, pstr0 = 0)
            }

            list(x = x, y = y)
}
}

