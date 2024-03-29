% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reflections.R
\name{sysabs}
\alias{sysabs}
\title{Locate systematic absences}
\usage{
sysabs(hkl, SG)
}
\arguments{
\item{hkl}{An mX3 matrix or a data frame whose rows are the
three integers corresponding to the Miller
indices.}

\item{SG}{A character. The extended Hermann-Mauguin symbol
of the crystallographic space group.}
}
\value{
idx A vector of integers corresponding to the
            position, in the array \code{mhkl}, in which the
            Miller indices ARE NOT systematically absent.
            The position of systematically-absent reflections
            can be found using !idx.
}
\description{
Given an mX3 matrix of Miller indices, this function returns
those indices corresponding to valid reflections, i.e. to
reflections which are not systematic absences.
}
\details{
Crystallography symmetry forces constraints on data from
x-ray diffraction. One of these constraints consists in the
full cancellation of reflections with certain Miller indices.
It is said that the reflection with that specific Miller index
is systematically absent. For example, in data corresponding
to a crystal with space group C 2, general reflections like
(h,k,l) must obey h+k=2n (even number). Thus, the Miller
indices (2,3,1) are a systematic absence because 2+3=5 (odd).
}
\examples{
# C 2 monoclinic space group (special setting)
csym <- cryst_symm(15,set=5)
print(csym$SG)

# Create a full data frame of Miller indices
hkl <- expand.grid(H=-4:4,K=-4:4,L=-4:4)

# Index corresponding to valid reflections
# (not systematic absences)
idx <- sysabs(hkl,csym$SG)

# Indices for all reflections
fulldx <- 1:length(hkl[,1])

# Index corresponding to systematic absences
jdx <- fulldx[-idx]

# A couple of systematic absences
hkl[jdx[1:2],]

}
