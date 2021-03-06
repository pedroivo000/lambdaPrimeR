% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thermodynamic_calculations.R
\name{nupack_pairs}
\alias{nupack_pairs}
\title{Calculate base-pairing probability using NUPACK's \code{pairs} tool.}
\usage{
nupack_pairs(input_sequence, melting_temperature, paired = TRUE,
  dangles = c("some", "none", "all"), cutoff = NULL)
}
\arguments{
\item{input_sequence}{A string corresponding to the input sequence.}

\item{melting_temperature}{A numeric value indicating the melting temperature
to be used during the NUPACK's calculations.}

\item{dangles}{A string with the way in which dangle energies are 
incorporated. Possible values: 

\itemize{
  \item \code{none}: No dangle energy are incorporated
  \item \code{some} (default): A dangle energy is incorporated for each 
  unpaired base flanking a duplex (a base flanking two duplexes contributes 
  only the minimum of the two possible dangle energies).
  \item \code{all}: A dangle energy is incorporated for each base flanking a 
  duplex regardless of whether it is paired.
}}

\item{cutoff}{A numeric value (default 0.001). Only probabilities above 
\code{cutoff} are written to output.}
}
\value{
A data frame with three columns: \code{p} corresponds to the 
probability of \code{ibase} being paired to \code{jbase}. If \code{paired} =
\code{FALSE}, column \code{p} contains the probability of 
\code{ibase} being unpaired.
}
\description{
Calculate base-pairing probability using NUPACK's \code{pairs} tool.
}
