\name{add_1}
\alias{add_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
add_1 to a binary number
}
\description{
j is a binary representation of a integer number.
The add_1 function adds 1 to the binary number j.
The input is j and the output is j+1 in binary representation.
}
\usage{
add_1(j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{a binary representation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
j+1    is returned in binary form
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hannah Lennon
}
\note{
This function requires the binary representation to be read in backwards!


E.g 1 in binary representation is 000001
    2 in binary representation is 000010
    3 in binary representation is 000011
    4 in binary representation is 000100
    5 in binary representation is 000101
    6 in binary representation is 000110
    However add_1(c(1,0,0,0,0,0)) = rev(000010)}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
c(1, 0, 1, 0, 0, 0) is 5 backwards so the result is 6 = 000110 backwards
}
\examples{
add_1(c(1, 0, 1, 0, 0, 0))


## The function is currently defined as
function (j)
{
    m <- length(j)
    if (!j[1])
        j[1] <- 1
    else {
        j[1] <- 0
        j[2:m] <- add_1(j[2:m])
    }
    return(j)
  }
}

