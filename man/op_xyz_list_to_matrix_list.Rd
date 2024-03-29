% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmetry.R
\name{op_xyz_list_to_matrix_list}
\alias{op_xyz_list_to_matrix_list}
\title{List of matrices and vectors of a specific space group}
\usage{
op_xyz_list_to_matrix_list(op_xyz_list)
}
\arguments{
\item{op_xyz_list}{A named list made of two vectors. The first vector, SYMOP, contains strings
describing the symmetry operators. The second vector, CENOP, contains strings describing the
centring of the unit cell.}
}
\value{
mat_ops_list A named list consisting of 3 lists. The first list, PG, contains
\eqn{3\times 3} point group matrices; the second list, T, contains the same number of
\eqn{3\times 1} translation vectors. The first matrix is always the identity matrix, the first
translation vector is always the null vector. The third list, C, consists of centering vectors;
the first centering vector is always the null vector. To summarize, the output looks like the
following:

           [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
where:
           I                = identity    3X3 matrix
           0                = null        3X1 vector
           M2,M3,...,Mn     = point group 3X3 matrices
           V2,V3,...,Cn     = translation 3X1 vectors
           C2,C3,...,Cm     = centering   3X1 vectors
}
\description{
Returns \eqn{3\times 3} matrices and \eqn{3\times 1} vectors corresponding to point group
operations, group translations and cell centring of a given space group.
}
\details{
A crystallographic space group consists of a series of transformations on a point
\eqn{(x_f,y_f,z_f)} in space that are mathematically implemented as the product of
a \eqn{3\times 3} point-group matrix and the point fractional coordinates, \eqn{(x_f,y_f,z_f)},
followed by a sum with a \eqn{3\times 1} translation vector. The complete set of points thus
produced can be cloned into a new and shifted set translated of an amount represented by a
\eqn{3\times 1} centring vector.
}
\examples{
# Symmetry operators for space group number 3, P 1 2 1
SG <- "P 1 2 1"
op_xyz_list <- syminfo_to_op_xyz_list(SG)
mat_ops_list <- op_xyz_list_to_matrix_list(op_xyz_list)
names(mat_ops_list)

}
