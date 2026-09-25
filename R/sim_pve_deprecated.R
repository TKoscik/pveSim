#' Execute Advanced 3D MRI Partial Volume Error (PVE) Simulation
#'
#' @description
#' Runs a comprehensive, multi-layered partial volume error simulation. The function maps 
#' continuous 3D geometries onto a discrete, customizable, anisotropic voxel grid, calculates 
#' mixed-tissue signal intensities based on structural and background contrast, layers on thermal 
#' MR noise, and returns both binary thresholded and linear mixel model volume measurements.
#' If requested, it outputs 3-plane cross-sectional slices through the grid center overlaid with 
#' a contrasting dotted contour indicating the analytical ground-truth shape boundary.
#'
#' @param shape A character string defining the target structure. Must be one of 
#'   \code{"sphere"}, \code{"ellipsoid"}, \code{"cylinder"}, \code{"ovoid"}, 
#'   \code{"irregular"}, or \code{"laminar"}.
#' @param params A named list containing the physical size parameters matching the chosen shape:
#'   \itemize{
#'     \item{\code{"sphere"}: \code{list(r)}
#'       \itemize{
#'         \item{\code{r}: Positive numeric value detailing the absolute radius of the sphere in voxel units.}
#'       }}
#'     \item{\code{"ellipsoid"}: \code{list(a, b, c)}
#'       \itemize{
#'         \item{\code{a}: Semi-major axis length along the X direction (in voxel units).}
#'         \item{\code{b}: Semi-minor axis length along the Y direction (in voxel units).}
#'         \item{\code{c}: Semi-minor axis length along the Z direction (in voxel units).}
#'       }}
#'     \item{\code{"cylinder"}: \code{list(r, h)}
#'       \itemize{
#'         \item{\code{r}: Cross-sectional radius of the cylinder base (in voxel units).}
#'         \item{\code{h}: Total longitudinal height/length of the cylinder barrel (in voxel units).}
#'       }}
#'     \item{\code{"ovoid"}: \code{list(a, b, c, k)}
#'       \itemize{
#'         \item{\code{a}: Base semi-axis scaling dimension along the X direction.}
#'         \item{\code{b}: Base semi-axis scaling dimension along the Y direction.}
#'         \item{\code{c}: Base semi-axis scaling dimension along the Z direction.}
#'         \item{\code{k}: Asymmetry scaling parameter. Values greater than 0 pinch one structural pole along the Z-axis while inflating the opposite pole to create a realistic egg/ovoid taper.}
#'       }}
#'     \item{\code{"irregular"}: \code{list(base_r, amp1, freq1, amp2, freq2)}
#'       \itemize{
#'         \item{\code{base_r}: Core underlying radius of the baseline unperturbed sphere.}
#'         \item{\code{amp1}: Physical amplitude (maximum height/depth in voxel units) of primary surface undulations.}
#'         \item{\code{freq1}: Spatial frequency (number of wave cycles per sphere sweep) of primary surface undulations.}
#'         \item{\code{amp2}: Physical amplitude of secondary overlapping surface textures.}
#'         \item{\code{freq2}: Spatial frequency of secondary overlapping surface textures.}
#'       }}
#'     \item{\code{"laminar"}: \code{list(base_z, curve_amp, curve_freq, thickness)}
#'       \itemize{
#'         \item{\code{base_z}: Absolute central vertical anchor plane elevation along the Z-axis.}
#'         \item{\code{curve_amp}: Maximum peak-to-trough amplitude height of folding sheet waves.}
#'         \item{\code{curve_freq}: Density or spatial frequency of folding sheet waves (higher values resemble highly gyrified cortical folds).}
#'         \item{\code{thickness}: Absolute continuous width/thickness of the tissue ribbon sheet layer (in voxel units).}
#'       }}
#'   }
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
#' @return A standard named list containing quantitative experimental metrics.
#' 
#' @export
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
  if (shape == "sphere") { max_dim <- params$r * 2
  } else if (shape == "ellipsoid") { max_dim <- max(params$a, params$b, params$c) * 2
  } else if (shape == "cylinder") { max_dim <- sqrt((2 * params$r)^2 + (params$h)^2) *2
  } else if (shape == "ovoid") { max_dim <- max(params$a, params$b, params$c) * 2
  } else if (shape == "irregular") { max_dim <- (params$base_r + params$amp1 + params$amp2) * 2
  } else if (shape == "laminar") { max_dim <- max(c(params$thickness, params$curve_amp * 2)) * 2
  } else { stop("Unknown shape designated.") }

  # Account for sub-voxel translation offsets shifting the shape out of center
  max_offset_shift <- max(abs(offset))
  # Compute total uniform cubic grid span required to hold the rotated/shifted shape
  total_needed_space <- max_dim + (2 * max_offset_shift) + buffer
  grid_span <- as.integer(ceiling(total_needed_space / min(res_inplane, res_z)))
  # Enforce an odd number so the voxel grid centers perfectly on (0,0,0)
  if (grid_span %% 2 == 0) { grid_span <- grid_span + 1 }
  grid <- make_voxel_grid(grid_span, res_inplane, res_inplane, res_z)
  grid_dim <- c(grid_span, grid_span, grid_span)
  
  #grid_span <- as.integer(ceiling((max_dim * 2 + buffer) / min(res_inplane, res_z)))
  #grid <- make_voxel_grid(grid_span, res_inplane, res_inplane, res_z)
  #grid_dim <- c(grid_span, grid_span, grid_span)
  
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
  #segmentation_threshold <- (int_structure + int_background) / 2
  #if (int_structure > int_background) {
  #  hard_segmentation <- voxel_intensity_observed >= segmentation_threshold
  #} else {
  #  hard_segmentation <- voxel_intensity_observed <= segmentation_threshold
  #}
  #v_hard <- sum(hard_segmentation) * grid$voxel_vol
  #err_hard <- ((v_hard - v_true) / v_true) * 100
  
  # Unbiased Hard Segmentation mapping to eliminate anisotropic threshold plateaus
  segmentation_threshold <- (int_structure + int_background) / 2
  # Calculate the raw distance from background normalized by the contrast gap
  scaled_intensity <- (voxel_intensity_observed - int_background) / (int_structure - int_background)
  # Create an unbiased hard segmentation mask 
  # (Values above 0.5 round up, values below round down)
  hard_segmentation <- scaled_intensity >= 0.5
  v_hard <- sum(hard_segmentation) * grid$voxel_vol
  err_hard <- ((v_hard - v_true) / v_true) * 100
  
  estimated_occupancy <- (voxel_intensity_observed - int_background) / (int_structure - int_background)
  estimated_occupancy[estimated_occupancy < 0] <- 0.0
  estimated_occupancy[estimated_occupancy > 1] <- 1.0
  v_mixel <- sum(estimated_occupancy) * grid$voxel_vol
  err_mixel <- ((v_mixel - v_true) / v_true) * 100
  
  boundary_voxels <- sum(voxel_occupancy > 0.0 & voxel_occupancy < 1.0)
  cnr <- abs(int_structure - int_background) / max(0.0001, noise_sd)
  
  # 7. Render 3-Plane Slice Imagery with Analytical Ground-Truth Overlay
  if (generate_plots) {
    mid_idx <- round(grid_span / 2)
    axial_slice    <- voxel_intensity_observed[, , mid_idx]
    coronal_slice  <- voxel_intensity_observed[, mid_idx, ]
    sagittal_slice <- voxel_intensity_observed[mid_idx, , ]
    
    gray_palette <- gray.colors(256, start = 0, end = 1)
    
    # Establish high-resolution evaluation grids specifically for clean contours
    contour_res <- 200
    x_dense <- seq(min(grid$x), max(grid$x), length.out = contour_res)
    y_dense <- seq(min(grid$y), max(grid$y), length.out = contour_res)
    z_dense <- seq(min(grid$z), max(grid$z), length.out = contour_res)
    
    # Calculate fixed center point slices based on current run parameters
    fixed_x <- grid$x[mid_idx]
    fixed_y <- grid$y[mid_idx]
    fixed_z <- grid$z[mid_idx]
    
    # 7a. Axial Analytical Contour Matrix
    axial_dense_grid <- expand.grid(X = x_dense, Y = y_dense)
    axial_pts <- matrix(c(axial_dense_grid$X - offset[1], 
                          axial_dense_grid$Y - offset[2], 
                          rep(fixed_z, nrow(axial_dense_grid)) - offset[3]), ncol = 3)
    axial_pts_mapped <- rotate_coords(axial_pts, -rot_rad[1], -rot_rad[2], -rot_rad[3])
    axial_contour_mask <- matrix(is_inside_shape(axial_pts_mapped, shape, params), 
                                 nrow = contour_res, ncol = contour_res)
    
    # 7b. Coronal Analytical Contour Matrix
    coronal_dense_grid <- expand.grid(X = x_dense, Z = z_dense)
    coronal_pts <- matrix(c(coronal_dense_grid$X - offset[1], 
                            rep(fixed_y, nrow(coronal_dense_grid)) - offset[2], 
                            coronal_dense_grid$Z - offset[3]), ncol = 3)
    coronal_pts_mapped <- rotate_coords(coronal_pts, -rot_rad[1], -rot_rad[2], -rot_rad[3])
    coronal_contour_mask <- matrix(is_inside_shape(coronal_pts_mapped, shape, params), 
                                   nrow = contour_res, ncol = contour_res)
    
    # 7c. Sagittal Analytical Contour Matrix
    sagittal_dense_grid <- expand.grid(Y = y_dense, Z = z_dense)
    sagittal_pts <- matrix(c(rep(fixed_x, nrow(sagittal_dense_grid)) - offset[1], 
                             sagittal_dense_grid$Y - offset[2], 
                             sagittal_dense_grid$Z - offset[3]), ncol = 3)
    sagittal_pts_mapped <- rotate_coords(sagittal_pts, -rot_rad[1], -rot_rad[2], -rot_rad[3])
    sagittal_contour_mask <- matrix(is_inside_shape(sagittal_pts_mapped, shape, params), 
                                    nrow = contour_res, ncol = contour_res)
    
    # Open File Device and Plot Images with Layered Contours
    absolute_intensity_range <- c(0, 255)
    png(filename = plot_filename, width = 1250, height = 420, res = 110)
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    
    # Plot Axial base view + Ground Truth contour
    image(grid$x, grid$y, axial_slice, col = gray_palette,
          zlim = absolute_intensity_range,
          main = paste("Axial (Z =", round(fixed_z, 1), ")"),
          xlab = "X", ylab = "Y",
          axes=FALSE, frame.plot = FALSE, asp = 1)
    contour(x_dense, y_dense, axial_contour_mask,
            levels = 0.5,col = "darkred", lty = 3, lwd = 2.5,
            add = TRUE, drawlabels = FALSE)
    
    # Plot Coronal base view + Ground Truth contour
    image(grid$x, grid$z, coronal_slice, col = gray_palette,
          zlim = absolute_intensity_range,
          main = paste("Coronal (Y =", round(fixed_y, 1), ")"),
          xlab = "X", ylab = "Z",
          axes=FALSE, frame.plot = FALSE, asp = 1)
    contour(x_dense, z_dense, coronal_contour_mask, levels = 0.5,
            col = "darkred", lty = 3, lwd = 2.5, add = TRUE, drawlabels = FALSE)
    
    # Plot Sagittal base view + Ground Truth contour
    image(grid$y, grid$z, sagittal_slice, col = gray_palette,
          zlim = absolute_intensity_range,
          main = paste("Sagittal (X =", round(fixed_x, 1), ")"),
          xlab = "Y", ylab = "Z",
          axes=FALSE, frame.plot = FALSE, asp = 1)
    contour(y_dense, z_dense, sagittal_contour_mask, levels = 0.5,
            col = "darkred", lty = 3, lwd = 2.5, add = TRUE, drawlabels = FALSE)
    
    dev.off()
  }
  # Return data block
  list(true_volume             = v_true,
       hard_volume             = v_hard,
       hard_error_pct          = err_hard,
       mixel_volume            = v_mixel,
       mixel_error_pct         = err_mixel,
       boundary_voxels_count   = boundary_voxels,
       contrast_to_noise_ratio = cnr,
       voxel_dimensions        = grid_dim)
}
