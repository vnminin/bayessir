\name{SMCwithModPossionTL}
\alias{SMCwithModPossionTL}
\title{Sequential Monte Carlo algorithm}
\usage{
SMCwithModPossionTL(observedCounts, observedTimes, parameters, alphaBreaks,
  startMeans, nparticles, popSize, reportingRate, deltatval, ncrit)
}
\arguments{
  \item{observedCounts}{Integer vector of observed counts}

  \item{observedTimes}{Numeric vector of observation times}

  \item{parameters}{Numeric vector of size three containing
  the current values of the infectious contact rate, the
  recovery rate, and the rate at which immunity is lost,
  respectively}

  \item{alphaBreaks}{List containing the time-varying
  environmental force of infection and the times at which
  the rate changes, broken up by simulation intervals}

  \item{startMeans}{Numeric vector of size two containing
  the means of initial distributions for the numbers of
  susceptible and infected individuals respectively}

  \item{nparticles}{Number of particles}

  \item{popSize}{Population size}

  \item{reportingRate}{Probability of infected individuals
  seeking treatment}

  \item{deltatval}{Initial value to use for tau in
  tau-leaping algorithm}

  \item{ncrit}{Critical size for modified tau-leaping
  algorithm; if the population of a compartment is lower
  than this number a single step algorithm is used until
  the population gets above the critical size.}
}
\value{
log-likelihood and estimate of last state
}
\description{
Sequential Monte Carlo algorithm
}

