\name{ModPoissonTL}
\alias{ModPoissonTL}
\title{Modified homogeneous Poisson tau-leaping algorithm}
\usage{
ModPoissonTL(startState, popSize, startTime, intervalLength, parameters, alpha,
  deltatint, ncrit)
}
\arguments{
  \item{startState}{Integer vector of size two containing
  the initial number of susceptible and infected
  individuals respectively}

  \item{popSize}{Population size}

  \item{startTime}{Start time of simulation}

  \item{intervalLength}{Length of simulation interval}

  \item{parameters}{parameters Numeric vector of size seven
  containing the current values of the infectious contact
  rate, the recovery rate, the rate at which immunity is
  lost, and the powers}

  \item{alpha}{The value of the time-varying environmental
  force of infection}

  \item{deltatint}{Initial value to use for tau in
  tau-leaping algorithm}

  \item{ncrit}{Critical number}
}
\value{
Vector containing number of susceptible and infected at the
end of the simulation interval
}
\description{
Modified Poisson tau-leaping algorithm. Uses a tau-leaping
algorithm to simulate the numbers of susceptible and
infected individuals forward in time when value of
environmental force of infection remains constant
}

