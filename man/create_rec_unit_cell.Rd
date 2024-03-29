% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cry_generics_S3.R
\name{create_rec_unit_cell}
\alias{create_rec_unit_cell}
\title{S3 generic to create rec_unit_cell objects}
\usage{
create_rec_unit_cell(ar, ...)
}
\arguments{
\item{ar}{An object or objects used to select a method. These can be
reciprocal unit cell parameters, an object of class
rec_unit_cell, etc.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
An object of class "rec_unit_cell". It is a named list of length 6 whose
        last three slots are of "angle" class.
}
\description{
The rec_unit_cell object can be created starting from specific objects, files, etc.
}
\examples{
# Create a rec_unit_cell in default (no arguments)
ruc <- create_rec_unit_cell()
print(ruc)
}
