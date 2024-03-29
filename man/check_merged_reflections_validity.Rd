% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validations.R
\name{check_merged_reflections_validity}
\alias{check_merged_reflections_validity}
\title{Validity and compatibility of a cry object of
class 'merged_reflections'}
\usage{
check_merged_reflections_validity(x, message = FALSE)
}
\arguments{
\item{x}{Object of class 'merged_reflections'.}

\item{message}{A logical variable. If TRUE, the function
prints a message on the errors, if any.}
}
\value{
ans A logical value. TRUE means that the input is
            a valid object of class'merged_reflections'.
}
\description{
An object of class 'merged_reflections' is a named list of
length 4:
\describe{
          \item{ruc}{An object of class "rec_unit_cell".}
          \item{csym}{An object of class "cryst_symm".}
          \item{records}{A data frame containing the data.}
          \item{dtypes}{A character vector containing the
                type of data (Miller indices, structure
                factors, etc).}
        }
Internal consistency must be displayed between the object
'ruc' and the object 'csym' because groups of crystallographic
symmetries are compatible only with certain unit cells (and,
accordingly, certain reciprocal cells).
It is not possible to check consistency between dtypes and
the nature of data in each column of the data frame 'records',
but a check about length of 'dtypes' and number of columns is
possible. Therefore, the user should pay attention to the
nature of his/her data. Also, merged reflection data, having
to be compatible with crystal symmetry, have to display the
appropriate systematic absences. Users interested in keeping
systematic absences in the object, might want to look at the
object of class "raw_reflections".
}
\examples{
# Create an object of class 'merged_reflections'
# (default ara data associated with a cubic cell)
x <- merged_reflections()

# Check its validity
check_merged_reflections_validity(x)

# Now change reciprocal unit cell (to triclinic)
uc <- unit_cell(10,20,30,30,50,70)
ruc <- create_rec_unit_cell(uc)
x$ruc <- ruc

# Check validity again
check_merged_reflections_validity(x)

}
