\name{simInit}
\alias{simInit}
\title{Sample from initial distributions}
\usage{
simInit(n, startState)
}
\arguments{
  \item{n}{Number of particles}

  \item{startState}{Numeric vector of size two containing
  the means of initial distributions for the numbers of
  susceptible and infected individuals respectively}
}
\value{
Matrix of dimension n by 2, where each row is the initial
state for a particle trajectory
}
\description{
Simulates the initial number of susceptible and infected
individuals from independent Poisson initial distributions.
}

