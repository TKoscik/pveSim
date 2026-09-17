\name{is_inside_shape}
\alias{is_inside_shape}
\title{Evaluate Point Intersection within 3D Analytical Geometries}
\usage{
  is_inside_shape(coords, shape, params)
}
\arguments{
  \item{coords}{A numeric matrix or data frame with exactly 3 columns representing
    the \eqn{(X, Y, Z)} coordinates of the transformation sampling points.}
  
  \item{shape}{A character string defining the target structure. Must be one of
    \code{"sphere"}, \code{"ellipsoid"}, \code{"cylinder"}, \code{"egg"},
    \code{"irregular"}, or \code{"laminar"}.}
  
  \item{params}{A named list containing the physical size parameters matching the chosen shape:
      \itemize{
        \item{\code{"sphere"}: \code{list(r)} where \code{r} is the radius.}
        \item{\code{"ellipsoid"}: \code{list(a, b, c)} defining the semi-axes lengths.}
        \item{\code{"cylinder"}: \code{list(r, h)} tracking cross-section radius and total height.}
        \item{\code{"ovoid"}: \code{list(a, b, c, k)} where \code{k} represents asymmetric scaling.}
        \item{\code{"irregular"}: \code{list(base_r, amp1, freq1, amp2, freq2)} for surface ripples.}
        \item{\code{"laminar"}: \code{list(base_z, curve_amp, curve_freq, thickness)} for layered tracking.}
      }}
}
\value{
  A logical vector of length \eqn{N} indicating \code{TRUE} for points located
  inside or on the boundary of the geometric structure, and \code{FALSE} otherwise.
}
\description{
  Evaluates whether an \eqn{N \times 3} matrix of continuous spatial points falls
  inside a specified target geometry. This function forms the core spatial lookup
  mechanism for voxelization, supporting standard geometries alongside organic
  and layered structures seen in clinical MRI morphometry.
}
\examples{
  # Evaluate if a coordinate point falls inside a sphere of radius 3
  sample_pts <- matrix(c(1, 1, 1, 5, 0, 0), ncol = 3, byrow = TRUE)
  sphere_check <- is_inside_shape(sample_pts, "sphere", list(r = 3.0))
  print(sphere_check) # Output: TRUE FALSE
}
