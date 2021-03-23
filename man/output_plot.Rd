\name{output_plot}
\alias{output_plot}
\title{Creates MCEM iteration plots
}
\description{Creates MCEM iteration plots pdf in current working directory}

\usage{
output_plot(output, initial)
}
\arguments{
  \item{output}{A data.frame of output, each column creates a new plot
}
  \item{initial}{ initial value given to the MCEM algorithm}
}
\details{
If initial is not specified, it defaults zero for each parameter}
\value{
pdf plot
}

\author{
Hannah Lennon
<drhannahlennon@gmail.com>
}
\note{
Save pdf of the plot in current working directory}

\examples{
data(ins.data)

results <- MCEM(ins.data[1:50], Monte_Carlo_samples=c(1, 5), p=1, q=1, initial=c(0,0), prec=0.01, stopping_criteria="MA_k", k=1, num_iterations=3)

output_plot(results$MCEM_Iterations)


}

