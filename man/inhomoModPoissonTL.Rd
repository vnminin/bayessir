\name{inhomoModPoissonTL}
\alias{inhomoModPoissonTL}
\title{Modified inhomogeneous Poisson tau-leaping algorithm}
\usage{
inhomoModPoissonTL(startState, popSize, startTime, intervalLength, parameters,
  alphas, allbreaks, deltatint, ncrit)
}
\arguments{
  \item{startState}{Integer vector of size two containing
  the initial number of susceptible and infected
  individuals respectively}

  \item{popSize}{Population size}

  \item{startTime}{Start time of simulation}

  \item{intervalLength}{Length of simulation interval}

  \item{parameters}{Numeric vector of size three containing
  the current values of the infectious contact rate, the
  recovery rate, and the rate at which immunity is lost,
  respectively}

  \item{alphas}{Numeric vector containing the values of
  alpha}

  \item{allbreaks}{Numeric vector containing the times at
  which the value of alpha changes}

  \item{deltatint}{Initial value to use for tau in
  tau-leaping algorithm}

  \item{ncrit}{Critical size for modified tau-leaping
  algorithm; if the population of a compartment is lower
  than this number a single step algorithm is used until
  the population gets above the critical size.}
}
\value{
Vector containing number of susceptible and infected at the
end of the simulation interval
}
\description{
Uses a tau-leaping algorithm to simulate the numbers of
susceptible and infected individuals forward in time when
value of environmental force of infection changes
}

