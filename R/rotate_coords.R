#' Apply 3D Rotation Matrix to Spatial Coordinates
#'
#' @description
#' Transforms an \eqn{N \times 3} matrix of continuous 3D spatial points using 
#' a combined Roll-Pitch-Yaw (\eqn{X-Y-Z}) rotation matrix. This function is 
#' primarily utilized by the simulation engine to apply random angular tilts to 
#' target geometries, mimicking arbitrary patient head positions inside an MRI scanner.
#'
#' @param coords A numeric matrix or data frame with exactly 3 columns representing 
#'   the \eqn{(X, Y, Z)} coordinates of points to transform. Row count \eqn{N} can be arbitrary.
#' @param rx A numeric value specifying the Roll angle around the X-axis in radians.
#' @param ry A numeric value specifying the Pitch angle around the Y-axis in radians.
#' @param rz A numeric value specifying the Yaw angle around the Z-axis in radians.
#'
#' @details
#' The rotation matrices are evaluated sequentially using a right-handed coordinate 
#' system. The unified rotation operator \eqn{R} is compiled via intrinsic matrix multiplication 
#' \eqn{R = R_z \times R_y \times R_x}. Coordinates are transformed by evaluating \eqn{Coords \times R^T}.
#' If all input angles are exactly \code{0}, the function skips multiplication and returns 
#' the baseline coordinate configuration instantly to preserve system resources.
#'
#' @return A numeric matrix of dimensions \eqn{N \times 3} containing the rotated spatial coordinates.
#' 
#' @export
#'
#' @examples
#' # Establish a sample coordinate space containing two distinct points
#' points <- matrix(c(1, 0, 0, 
#'                     0, 1, 0), ncol = 3, byrow = TRUE)
#' 
#' # Apply a 90-degree (pi/2 radians) yaw rotation around the Z-axis
#' rotated_points <- rotate_coords(points, rx = 0, ry = 0, rz = pi/2)
#' print(rotated_points)
rotate_coords <- function(coords, rx, ry, rz) {
  
  # Ensure spatial matrix properties match dimensions explicitly
  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop("Argument 'coords' must be a structured numeric matrix or data frame.")
  }
  if (ncol(coords) != 3) {
    stop("Argument 'coords' must possess exactly 3 columns representing X, Y, and Z axes.")
  }
  if (nrow(coords) == 0) {
    return(matrix(numeric(0), ncol = 3))
  }
  
  # Convert data frames to strict numeric matrices for matrix multiplication speeds
  if (is.data.frame(coords)) {
    coords <- as.matrix(coords)
  }
  
  # Strip missing entries to guard linear arithmetic combinations
  if (any(is.na(coords)) || any(is.na(c(rx, ry, rz)))) {
    stop("Input matrices and rotation angles cannot contain missing values (NA).")
  }
  
  # Resource shortcut: If no rotation is demanded, bypass system operators
  if (rx == 0 && ry == 0 && rz == 0) {
    return(coords)
  }
  
  # Compute direction operators explicitly
  cx <- cos(rx); sx <- sin(rx)
  cy <- cos(ry); sy <- sin(ry)
  cz <- cos(rz); sz <- sin(rz)
  
  # Build individual 3D Euler element transformation parameters
  Rx <- matrix(c(1,   0,   0,
                 0,  cx, -sx,
                 0,  sx,  cx), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Pitch adjustments along the secondary axis tracking configuration
  Ry <- matrix(c( cy,  0,  sy,
                  0,  1,   0,
                  -sy,  0,  cy), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Yaw adjustments capturing orientation sweeping planes
  Rz <- matrix(c(cz, -sz,   0,
                 sz,  cz,   0,
                 0,   0,   1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Compile rotation matrices and apply transposed coordinate mappings
  R <- Rz %*% (Ry %*% Rx)
  
  return(coords %*% t(R))
}
