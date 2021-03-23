\name{binary}
\alias{binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{binary
}
\description{binary
Gives the binary representaion of value q in a vector of lenghth n }
\usage{
binary(q, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
a numeric value
}
  \item{n}{
a numeric value}
}
\details{
If n is not specified, it defaults to a output a vector of length 10, so q must be less than 1024}
\value{
a vector of length n in binary representation
}

\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}


\examples{
binary(10)
binary(10, 5)
}
