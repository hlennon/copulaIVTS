\name{print_iteration}
\alias{print_iteration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{print_iteration
}
\description{print_iteration
Gives the updated parameter values in correct form }
\usage{
print_iteration(para, p, q, diff, finish, t, m)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{para}{vector of updated parameter values }
  \item{p}{ar order  }
  \item{q}{ma order  }
  \item{diff}{difference from last iteration }
    \item{method}{"Delta_l_x" or "MA_1" or "MA_k" : Stopping Criteria used }
  \item{t}{index of iteration  }
  \item{m}{number of Monte-Carlo samples used in eevaluating the E-step}
}
\details{
None

}

\value{
prints updated values of the iterations
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
print_iteration(c(0.4, 0.4, 0.4), 2, 2, 4, TRUE, 10,  1)



## The function is currently defined as
print_iteration<-function(para, p, q, diff, finish, t, m){


            if(p==1){
                        cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Delta_x=',
                            diff, 'Stop?=', finish,  "m=", m,'\n')
            }

            if(p==2){
                        cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Delta_x=',
                            diff, 'Stop?=', finish,  "m=", m,'\n')
            }


}

}

