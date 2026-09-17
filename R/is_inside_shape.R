#' Evaluate Point Intersection within 3D Analytical Geometries
#'
#' @description
#' Evaluates whether an \eqn{N \times 3} matrix of continuous spatial points falls 
#' inside a specified target geometry. This function forms the core spatial lookup 
#' mechanism for voxelization, supporting standard geometries alongside organic 
#' and layered structures seen in clinical MRI morphometry.
#'
#' @param coords A numeric matrix or data frame with exactly 3 columns representing 
#'   the \eqn{(X, Y, Z)} coordinates of the transformation sampling points.
#' @param shape A character string defining the target structure. Must be one of 
#'   \code{"sphere"}, \code{"ellipsoid"}, \code{"cylinder"}, \code{"egg"}, 
#'   \code{"irregular"}, or \code{"laminar"}.
#' @param params A named list containing the physical size parameters matching the chosen shape:
#'   \itemize{
#'     \item{\code{"sphere"}: \code{list(r)} where \code{r} is the radius.}
#'     \item{\code{"ellipsoid"}: \code{list(a, b, c)} defining the semi-axes lengths.}
#'     \item{\code{"cylinder"}: \code{list(r, h)} tracking cross-section radius and total height.}
#'     \item{\code{"ovoid"}: \code{list(a, b, c, k)} where \code{k} represents asymmetric scaling.}
#'     \item{\code{"irregular"}: \code{list(base_r, amp1, freq1, amp2, freq2)} for surface ripples.}
#'     \item{\code{"laminar"}: \code{list(base_z, curve_amp, curve_freq, thickness)} for layered tracking.}
#'   }
#'
#' @return A logical vector of length \eqn{N} indicating \code{TRUE} for points located 
#'   inside or on the boundary of the geometric structure, and \code{FALSE} otherwise.
#' 
#' @export
#'
#' @examples
#' # Evaluate if a coordinate point falls inside a sphere of radius 3
#' sample_pts <- matrix(c(1, 1, 1, 5, 0, 0), ncol = 3, byrow = TRUE)
#' sphere_check <- is_inside_shape(sample_pts, "sphere", list(r = 3.0))
#' print(sphere_check) # Output: TRUE FALSE
is_inside_shape <- function(coords, shape, params) {
  
  # 1. Structural Checks and Coercions
  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop("Argument 'coords' must be a structured numeric matrix or data frame.")
  }
  if (ncol(coords) != 3) {
    stop("Argument 'coords' must possess exactly 3 columns representing X, Y, and Z axes.")
  }
  if (is.data.frame(coords)) {
    coords <- as.matrix(coords)
  }
  if (nrow(coords) == 0) {
    return(logical(0))
  }
  
  valid_shapes <- c("sphere", "ellipsoid", "cylinder", "ovoid", "irregular", "laminar")
  if (!(shape %in% valid_shapes)) {
    stop(paste("Shape profile must be one of:", paste(valid_shapes, collapse = ", ")))
  }
  if (!is.list(params)) {
    stop("Argument 'params' must be a named list containing structural dimension parameters.")
  }
  
  # Deconstruct components for vectorized processing
  x <- coords[, 1]
  y <- coords[, 2]
  z <- coords[, 3]
  
  # 2. Shape Evaluations
  if (shape == "sphere") {
    if (is.null(params$r) || params$r <= 0) stop("Sphere parameter 'r' must be positive.")
    return((x^2 + y^2 + z^2) <= params$r^2)
    
  } else if (shape == "ellipsoid") {
    if (any(c(is.null(params$a), is.null(params$b), is.null(params$c))) || 
        any(c(params$a, params$b, params$c) <= 0)) {
      stop("Ellipsoid parameters 'a', 'b', and 'c' must be positive.")
    }
    return((x / params$a)^2 + (y / params$b)^2 + (z / params$c)^2 <= 1)
    
  } else if (shape == "cylinder") {
    if (any(c(is.null(params$r), is.null(params$h))) || any(c(params$r, params$h) <= 0)) {
      stop("Cylinder parameters 'r' and 'h' must be positive.")
    }
    return((x^2 + y^2 <= params$r^2) & (abs(z) <= (params$h / 2)))
    
  } else if (shape == "ovoid") {
    if (any(c(is.null(params$a), is.null(params$b), is.null(params$c), is.null(params$k))) ||
        any(c(params$a, params$b, params$c) <= 0)) {
      stop("Egg parameters 'a', 'b', and 'c' must be positive numeric values.")
    }
    # Calculate asymmetric ovoid scaling profile along the axis
    scaled_radius_sq <- (1 + params$k * z)
    scaled_radius_sq[scaled_radius_sq <= 0] <- 0.0001 # Edge distortion limit
    return(((x / params$a)^2 + (y / params$b)^2) / scaled_radius_sq + (z / params$c)^2 <= 1)
    
  } else if (shape == "irregular") {
    # Check for spherical harmonic perturbation parameters
    req_fields <- c("base_r", "amp1", "freq1", "amp2", "freq2")
    if (any(sapply(params[req_fields], is.null))) {
      stop("Irregular shape list requires: base_r, amp1, freq1, amp2, and freq2 parameters.")
    }
    
    r_coords <- sqrt(x^2 + y^2 + z^2)
    r_coords[r_coords == 0] <- 0.0001 # Safe division at origin coordinates
    
    theta <- acos(z / r_coords)
    phi <- atan2(y, x)
    
    r_boundary <- params$base_r + 
      params$amp1 * sin(params$freq1 * theta) * cos(params$freq1 * phi) +
      params$amp2 * cos(params$freq2 * theta)
    
    return(r_coords <= r_boundary)
    
  } else if (shape == "laminar") {
    req_fields <- c("base_z", "curve_amp", "curve_freq", "thickness")
    if (any(sapply(params[req_fields], is.null)) || params$thickness <= 0) {
      stop("Laminar shape list requires positive 'thickness' and numerical curve bounds.")
    }
    
    # Evaluate localized sheet centerline profile height
    z_profile <- params$base_z + params$curve_amp * sin(params$curve_freq * x) * cos(params$curve_freq * y)
    return(abs(z - z_profile) <= (params$thickness / 2))
  }
}
