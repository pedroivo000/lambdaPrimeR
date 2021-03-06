% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primer_design.R
\name{get_overlaps}
\alias{get_overlaps}
\title{Extract overlaps from input sequences.}
\usage{
get_overlaps(sequence_inputs, position = NULL, length = NULL)
}
\arguments{
\item{sequence_inputs}{A dataframe containing the target and template input
sequences, generated by \code{read_sequences}.}

\item{position}{Numeric value indicating the insertion coordinate}

\item{length}{Length of the overlaping region to be extract from the input
sequences. Defaults to different values depending on the value of \code{type}:

\itemize{
  \item \code{type} = 'template' => \code{length} = 20
  \item \code{type} = 'target' => \code{length} = 15
}}
}
\value{
A dataframe containing the overlap sequence and orientation for each
input type.
}
\description{
Extract overlaps from input sequences.
}
\keyword{internal}
