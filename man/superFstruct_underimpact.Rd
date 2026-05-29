\name{superFstruct_underimpact}
\alias{superFstruct_underimpact}
\docType{data}
\title{Binary supervision structure to reconstruct the issue of underimpact of partial supervision.}
\description{
  This dataset provides a concrete superivison structure:
  - `superF` matrix of size 150 x 3 with partial supervision,
  - `ind` vector with indices of unsupervised observations,
  - `tind` vector with indicies of observations selected to be in the test dataset,
  - `tclass` vector with class membership of the observations selected
  to be in the test dataset.
  
  This supervision structure is meant to reproduce a particular realization of
  phenomenon of underimpact of partial supervision specific to the iris dataset.
}
\usage{data(superFstruct_underimpact)}
\format{A list with: a matrix of size 150 x 3, and three vectors.}
\keyword{datasets}