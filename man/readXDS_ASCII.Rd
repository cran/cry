% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readXDS_ASCII.R
\name{readXDS_ASCII}
\alias{readXDS_ASCII}
\title{Load an XDS_ASCII file.}
\usage{
readXDS_ASCII(filename, message = FALSE)
}
\arguments{
\item{filename}{A character string. The path to a valid
XDS ASCII file.}

\item{message}{A logical variable. If TRUE (default) the
function prints a message highlighting what is
included in the xds header. If \code{filename} is not a valid XDS
ascii file, the function returns `NULL` and
prints out a warning message.}
}
\value{
A named list (see details).
}
\description{
Function to load XDS_ASCII.HKL files into a named list with
three components called \emph{processing_info}, \emph{header}
and \emph{reflections} (see details further down).
}
\details{
This function reads in all data from an XDS_ASCII data file
and organises them into a named list. The list's name are:
\describe{
  \item{processing_info}{This list component includes three
        logical variables, MERGE, FRIEDEL and PROFILE. Their
        TRUE/FALSE value reflect features of the XDS_ASCII
        file connected with the specific processing performed
        to obtain the file itself
        (for more details see \url{https://xds.mr.mpg.de/}).}
   \item{header}{This list includes several components, like
                 for instance SPACE_GROUP_NUMBER or
                 UNIT_CELL_CONSTANTS, which give informations
                 on the crystal and the experiment generating
                 the data.}
   \item{reflections}{This data.frame includes the actual
                     experimental data, i.e. the observations
                     collected during the X-ray diffraction
                     experiment on the crystal (or crystals).
                     The number and type of columns can vary.}
}
}
\examples{

# Load one of the XDS ASCII files included with
# this distribution of cry
datadir <- system.file("extdata",package="cry")
filename <- file.path(datadir,"xds00_ascii.hkl")
ltmp <- readXDS_ASCII(filename,message=FALSE)
print(names(ltmp))
print(ltmp$reflections[1:5,])

}
