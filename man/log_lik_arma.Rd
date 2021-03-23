\name{log_lik_arma}
\alias{log_lik_arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-likelihood function of continuous Gaussian ARMA  model}
\description{
Direct estimation of the log-likelihood function of continuous Gaussian ARMA model}
\usage{
log_lik_arma(xxx, para, order)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxx}{
%%     ~~Describe \code{xxx} here~~
}
  \item{para}{
%%     ~~Describe \code{para} here~~
}
  \item{order}{
%%     ~~Describe \code{order} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## The function is currently defined as
function (xxx, para, order)
{
    m <- ncol(xxx)
    Ij <- numeric(m)
    for (i in 1:m) {
        Ij[i] <- exp(arma_loglikelihood(para, xxx[, i], ord = order))
    }
    return(Ij)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
