\name{rotate_coords}
\alias{rotate_coords}
\title{Apply 3D Rotation Matrix to Spatial Coordinates}
\usage{
  rotate_coords(coords, rx, ry, rz)
}
\arguments{
  \item{coords}{A numeric matrix or data frame with exactly 3 columns representing
    the \eqn{(X, Y, Z)} coordinates of points to transform. Row count \eqn{N} can be arbitrary.}
  
  \item{rx}{A numeric value specifying the Roll angle around the X-axis in radians.}
  
  \item{ry}{A numeric value specifying the Pitch angle around the Y-axis in radians.}
  
  \item{rz}{A numeric value specifying the Yaw angle around the Z-axis in radians.}
}
\value{
  A numeric matrix of dimensions \eqn{N \times 3} containing the rotated spatial coordinates.
}
\description{
  Transforms an \eqn{N \times 3} matrix of continuous 3D spatial points using
  a combined Roll-Pitch-Yaw (\eqn{X-Y-Z}) rotation matrix. This function is
  primarily utilized by the simulation engine to apply random angular tilts to
  target geometries, mimicking arbitrary patient head positions inside an MRI scanner.
}
\details{
  The rotation matrices are evaluated sequentially using a right-handed coordinate
  system. The unified rotation operator \eqn{R} is compiled via intrinsic matrix multiplication
  \eqn{R = R_z \times R_y \times R_x}. Coordinates are transformed by evaluating \eqn{Coords \times R^T}.
  If all input angles are exactly \code{0}, the function skips multiplication and returns
  the baseline coordinate configuration instantly to preserve system resources.
}
\examples{
  # Establish a sample coordinate space containing two distinct points
  points <- matrix(c(1, 0, 0,
                     0, 1, 0), ncol = 3, byrow = TRUE)
  
  # Apply a 90-degree (pi/2 radians) yaw rotation around the Z-axis
  rotated_points <- rotate_coords(points, rx = 0, ry = 0, rz = pi/2)
  print(rotated_points)
}
