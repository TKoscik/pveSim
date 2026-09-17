#' Execute Advanced 3D MRI Partial Volume Error (PVE) Simulation
#'
#' @description
#' Runs a comprehensive, multi-layered partial volume error simulation. The function maps 
#' continuous 3D geometries onto a discrete, customizable, anisotropic voxel grid, calculates 
#' mixed-tissue signal intensities based on structural and background contrast, layers on thermal 
#' MR noise, and returns both binary thresholded and linear mixel model volume measurements.
#'
#' @param shape A character string defining the target structure. Must be one of 
#'   \code{"sphere"}, \code{"ellipsoid"}, \code{"cylinder"}, \code{"egg"}, 
#'   \code{"irregular"}, or \code{"laminar"}.
#' @param params A named list containing the physical size parameters matching the chosen shape.
#' @param res_inplane A positive numeric value detailing the in-plane X and Y spatial resolution. 
#'   Defaults to \code{1.0}.
#' @param res_z A positive numeric value detailing the slice thickness (Z resolution). 
#'   Defaults to \code{3.0}.
#' @param offset A numeric vector of length 3 detailing the sub-voxel \code{c(x, y, z)} 
#'   translation displacement offsets. Defaults to \code{c(0, 0, 0)}.
#' @param rot_deg A numeric vector of length 3 detailing the roll, pitch, and yaw 
#'   rotation angles in degrees. Defaults to \code{c(0, 0, 0)}.
#' @param int_structure A numeric value defining the pure baseline MR signal intensity 
#'   of the structure. Defaults to \code{100}.
#' @param int_background A numeric value defining the pure baseline MR signal intensity 
#'   of the background tissue. Defaults to \code{50}.
#' @param noise_sd A non-negative numeric value defining the standard deviation of 
#'   the ambient Gaussian electronic/thermal MR noise. Defaults to \code{2.0}.
#' @param subsamples An integer value controlling the supersampling resolution (subdivisions per voxel axis) 
#'   used for numerical fractional integration. Total sub-voxels evaluated per voxel is \code{subsamples^3}. 
#'   Defaults to \code{4}.
#' @param generate_plots A logical flag to manually trigger PNG slice generation through the spatial 
#'   center of the voxel matrix array. Defaults to \code{FALSE}.
#' @param plot_filename A character string specifying the file pathway where the diagnostic 
#'   cross-section PNG layout should be exported if \code{generate_plots = TRUE}. 
#'   Defaults to \code{"mri_slice_output.png"}.
#'
#' @return A standard named list containing the following quantitative experimental metrics:
#' \describe{
#'   \item{\code{true_volume}}{The calculated baseline ground-truth continuous volume footprint.}
#'   \item{\code{hard_volume}}{The estimated volume calculated via traditional binary 50\% threshold segmentation.}
#'   \item{\code{hard_error_pct}}{The resulting percentage error under hard thresholding (positive values indicate overestimation, negative indicate underestimation).}
#'   \item{\code{mixel_volume}}{The estimated volume calculated using soft fractional un-mixing (The Linear Mixel Model).}
#'   \item{\code{mixel_error_pct}}{The resulting percentage error under soft fractional un-mixing.}
#'   \item{\code{boundary_voxels_count}}{The absolute number of voxels holding a blended tissue profile.}
#'   \item{\code{contrast_to_noise_ratio}}{The explicit Contrast-to-Noise Ratio (CNR) profile of the simulation run.}
#'   \item{\code{voxel_dimensions}}{A vector tracking the dimensions of the generated simulation grid array.}
#' }
#' 
#' @export
#'
#' @examples
#' # Simulate an anisotropic acquisition of an ellipsoid structure with a 15-degree tilt
#' my_sim <- sim_pve(
#'   shape = "ellipsoid", 
#'   params = list(a = 5.0, b = 3.5, c = 2.5),
#'   res_inplane = 1.0, res_z = 2.5,
#'   rot_deg = c(15, 0, 0),
#'   generate_plots = FALSE
#' )
#' print(my_sim$hard_error_pct)
sim_pve <- function(shape = "sphere", params = list(r = 4.0), 
                    res_inplane = 1.0, res_z = 3.0, 
                    offset = c(0.0, 0.0, 0.0), rot_deg = c(0.0, 0.0, 0.0),
                    int_structure = 100, int_background = 50, noise_sd = 2.0,
                    subsamples = 4, 
                    generate_plots = FALSE, plot_filename = "mri_slice_output.png") {
  
  # 1. Parameter Boundary Verifications
  if (res_inplane <= 0 || res_z <= 0) stop("Grid resolutions must be positive numbers.")
  if (length(offset) != 3 || length(rot_deg) != 3) stop("Offset and rot_deg must be numeric vectors of length 3.")
  if (noise_sd < 0) stop("Noise standard deviation cannot be a negative value.")
  if (int_structure == int_background) stop("Structure intensity and background intensity cannot be identical.")
  if (subsamples <= 0 || subsamples %% 1 != 0) stop("Subsamples parameter must be a strictly positive integer.")
  
  # 2. Derive Dynamic Bounding Box Dimensions
  buffer <- max(c(res_inplane, res_z)) * 4
  if (shape == "sphere") { max_dim <- params$r
  } else if (shape == "ellipsoid") { max_dim <- max(params$a, params$b, params$c)
  } else if (shape == "cylinder") { max_dim <- max(params$r, params$h/2)
  } else if (shape == "egg") { max_dim <- max(params$a, params$b, params$c)
  } else if (shape == "irregular") { max_dim <- params$base_r + params$amp1 + params$amp2
  } else if (shape == "laminar") { max_dim <- max(c(params$thickness, params$curve_amp * 2)) 
  } else { stop("Unknown shape designated.") }
  
  grid_span <- as.integer(ceiling((max_dim * 2 + buffer) / min(res_inplane, res_z)))
  grid <- make_voxel_grid(grid_span, res_inplane, res_inplane, res_z)
  grid_dim <- c(grid_span, grid_span, grid_span)
  
  # 3. Setup Numeric Supersampling Map
  rot_rad <- rot_deg * pi / 180
  sub_offsets <- seq(-0.5 + 1/(2*subsamples), 0.5 - 1/(2*subsamples), length.out = subsamples)
  sub_grid <- expand.grid(dx = sub_offsets * res_inplane, 
                          dy = sub_offsets * res_inplane, 
                          dz = sub_offsets * res_z)
  n_sub <- nrow(sub_grid)
  sub_grid_matrix <- as.matrix(sub_grid)
  
  voxel_occupancy <- array(0.0, dim = grid_dim)
  
  # 4. Core Voxelization Loop
  for (i in seq_len(grid_span)) {
    for (j in seq_len(grid_span)) {
      for (k in seq_len(grid_span)) {
        v_center <- c(grid$x[i], grid$y[j], grid$z[k])
        
        # Build evaluation blocks for continuous matrix checks
        pts <- matrix(rep(v_center, each = n_sub), ncol = 3) + sub_grid_matrix
        pts_transformed <- matrix(c(pts[,1] - offset[1], 
                                    pts[,2] - offset[2], 
                                    pts[,3] - offset[3]), ncol = 3)
        
        pts_transformed <- rotate_coords(pts_transformed, -rot_rad[1], -rot_rad[2], -rot_rad[3])
        inside_flags <- is_inside_shape(pts_transformed, shape, params)
        voxel_occupancy[i, j, k] <- sum(inside_flags) / n_sub
      }
    }
  }
  
  v_true <- sum(voxel_occupancy) * grid$voxel_vol
  if (v_true == 0) v_true <- 0.0001
  
  # 5. Contrast and Instrument Noise Modeling
  voxel_intensity_pure <- voxel_occupancy * int_structure + (1.0 - voxel_occupancy) * int_background
  noise_matrix <- array(rnorm(prod(grid_dim), mean = 0, sd = noise_sd), dim = grid_dim)
  voxel_intensity_observed <- voxel_intensity_pure + noise_matrix
  
  # 6. Extraction Segmentations
  segmentation_threshold <- (int_structure + int_background) / 2
  if (int_structure > int_background) {
    hard_segmentation <- voxel_intensity_observed >= segmentation_threshold
  } else {
    hard_segmentation <- voxel_intensity_observed <= segmentation_threshold
  }
  v_hard <- sum(hard_segmentation) * grid$voxel_vol
  err_hard <- ((v_hard - v_true) / v_true) * 100
  
  estimated_occupancy <- (voxel_intensity_observed - int_background) / (int_structure - int_background)
  estimated_occupancy[estimated_occupancy < 0] <- 0.0
  estimated_occupancy[estimated_occupancy > 1] <- 1.0
  v_mixel <- sum(estimated_occupancy) * grid$voxel_vol
  err_mixel <- ((v_mixel - v_true) / v_true) * 100
  
  boundary_voxels <- sum(voxel_occupancy > 0.0 & voxel_occupancy < 1.0)
  cnr <- abs(int_structure - int_background) / max(0.0001, noise_sd)
  
  # 7. Render 3-Plane Slice Imagery (PNG)
  if (generate_plots) {
    mid_idx <- round(grid_span / 2)
    axial_slice    <- voxel_intensity_observed[, , mid_idx]
    coronal_slice  <- voxel_intensity_observed[, mid_idx, ]
    sagittal_slice <- voxel_intensity_observed[mid_idx, , ]
    
    gray_palette <- gray.colors(256, start = 0, end = 1)
    
    png(filename = plot_filename, width = 1200, height = 400, res = 100)
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    
    image(grid$x, grid$y, axial_slice, col = gray_palette, 
          main = paste("Axial Slice (Z =", round(grid$z[mid_idx], 1), ")"),
          xlab = "X Dimension", ylab = "Y Dimension", asp = 1)
    
    image(grid$x, grid$z, coronal_slice, col = gray_palette, 
          main = paste("Coronal Slice (Y =", round(grid$y[mid_idx], 1), ")"),
          xlab = "X Dimension", ylab = "Z Slice Axis", asp = 1)
    
    image(grid$y, grid$z, sagittal_slice, col = gray_palette, 
          main = paste("Sagittal Slice (X =", round(grid$x[mid_idx], 1), ")"),
          xlab = "Y Dimension", ylab = "Z Slice Axis", asp = 1)
    
    dev.off()
  }
  
  # Return data block
  list(
    true_volume             = v_true,
    hard_volume             = v_hard,
    hard_error_pct          = err_hard,
    mixel_volume            = v_mixel,
    mixel_error_pct         = err_mixel,
    boundary_voxels_count   = boundary_voxels,
    contrast_to_noise_ratio = cnr,
    voxel_dimensions        = grid_dim
  )
}
