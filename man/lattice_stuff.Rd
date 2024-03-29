% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clattice.R
\name{lattice_stuff}
\alias{lattice_stuff}
\title{Calculation of useful lattice parameters}
\usage{
lattice_stuff(a, b, c, aa, bb, cc)
}
\arguments{
\item{a}{A real number. One of the unit cell's side lengths, in angstroms.}

\item{b}{A real number. One of the unit cell's side lengths, in angstroms.}

\item{c}{A real number. One of the unit cell's side lengths, in angstroms.}

\item{aa}{A real number. One of the unit cell's angles, in degrees.}

\item{bb}{A real number. One of the unit cell's angles, in degrees.}

\item{cc}{A real number. One of the unit cell's angles, in degrees.}
}
\value{
A named vector of real numbers and length 16. The names are:
 \itemize{
  \item{sa.   Sine(aa)}
  \item{sb.   Sine(bb)}
  \item{sc.   Sine(cc)}
  \item{ca.   Cosine(aa)}
  \item{cb.   Cosine(bb)}
  \item{cc.   Cosine(cc)}
  \item{ar.   a* (reciprocal cell side length)}
  \item{br.   b* (reciprocal cell side length)}
  \item{cr.   c* (reciprocal cell side length)}
  \item{sar.  Sine of a reciprocal cell angle}
  \item{sbr.  Sine of a reciprocal cell angle}
  \item{scr.  Sine of a reciprocal cell angle}
  \item{car.  Cosine of a reciprocal cell angle}
  \item{cbr.  Cosine of a reciprocal cell angle}
  \item{ccr.  Cosine of a reciprocal cell angle}
  \item{V.    Volume of the unit cell in cubic angstroms}
 }
}
\description{
Calculation of useful lattice parameters
}
\examples{
datadir <- system.file("extdata",package="cry")
fname <- file.path(datadir,"1dei_phases.mtz")
hdr <- readMTZHeader(fname,message=FALSE)
ucell <- hdr$CELL
vtmp <- lattice_stuff(ucell[1],ucell[2],ucell[3],ucell[4],ucell[5],ucell[6])
vtmp[1:3]
vtmp[4:6]
vtmp[7:9]
vtmp[10:12]
vtmp[13:15]
vtmp[16]
}
