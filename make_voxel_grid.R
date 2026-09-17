#' Create a 3D Voxel Coordinate Grid System
#'
#' @description
#' Generates an structured, zero-centered 3D grid array containing coordinates 
#' for every individual voxel midpoint. This function serves as the underlying 
#' spatial reference grid for evaluating geometric object intersections and partial 
#' volume calculations under varying levels of voxel anisotropy.
#'
#' @param grid_size An integer value specifying the number of voxels along each 
#'   dimension of the cubic grid. The total number of voxels evaluated will be \code{grid_size^3}.
#' @param res_x A positive numeric value detailing the in-plane X resolution (voxel width). 
#'   Defaults to \code{1.0}.
#' @param res_y A positive numeric value detailing the in-plane Y resolution (voxel height). 
#'   Defaults to \code{1.0}.
#' @param res_z A positive numeric value detailing the slice thickness (Z-axis resolution). 
#'   Defaults to \code{3.0} to reflect realistic clinical MRI slice anisotropy.
#'
#' @return A standard named list containing the following geometric elements:
#' \describe{
#'   \item{\code{x}}{A numeric vector containing the midpoint coordinates along the X axis.}
#'   \item{\code{y}}{A numeric vector containing the midpoint coordinates along the Y axis.}
#'   \item{\code{z}}{A numeric vector containing the midpoint coordinates along the Z axis.}
#'   \item{\code{voxel_vol}}{The calculated single-voxel physical volume (\code{res_x * res_y * res_z}).}
#'   \item{\code{res}}{A numeric vector of length 3 tracking the spatial resolutions used (\code{c(res_x, res_y, res_z)}).}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate a standard clinical anisotropic MRI grid profile
#' my_grid <- make_voxel_grid(grid_size = 21, res_x = 1.0, res_y = 1.0, res_z = 3.0)
#' 
#' # Extract total single-voxel volume capacity
#' print(my_grid$voxel_vol)

make_voxel_grid <- function(grid_size, res_x = 1.0, res_y = 1.0, res_z = 3.0) {
  
  # Structural Input Validations
  if (grid_size <= 0 || grid_size %% 1 != 0) {
    stop("Argument 'grid_size' must be a strictly positive integer value.")
  }
  if (res_x <= 0 || res_y <= 0 || res_z <= 0) {
    stop("Voxel spatial resolutions (res_x, res_y, res_z) must be positive values.")
  }
  
  # Calculate symmetrical grid midpoint sequences centered precisely around (0,0,0)
  x_seq <- (seq_len(grid_size) - (grid_size + 1) / 2) * res_x
  y_seq <- (seq_len(grid_size) - (grid_size + 1) / 2) * res_y
  z_seq <- (seq_len(grid_size) - (grid_size + 1) / 2) * res_z
  
  # Compile structured geometry output list
  list(
    x         = x_seq,
    y         = y_seq,
    z         = z_seq, 
    voxel_vol = res_x * res_y * res_z,
    res       = c(res_x, res_y, res_z)
  )
}
