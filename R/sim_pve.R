#' Execute Advanced 3D MRI Partial Volume Error (PVE) Simulation
#'
#' @param shape Character string: "sphere", "ellipsoid", "cylinder", "ovoid", "irregular", or "laminar".
#' @param params Named list containing shape parameters (e.g. list(r = 4.0)).
#' @param res_vox Numeric vector of length 3 detailing physical voxel dimensions c(dx, dy, dz).
#' @param offset Sub-voxel translation offset vector c(x, y, z).
#' @param rot_deg Rotation angles in degrees c(roll, pitch, yaw).
#' @param int_structure Baseline MR signal intensity of structure.
#' @param int_background Baseline MR signal intensity of background.
#' @param noise_sd Standard deviation of Gaussian noise.
#' @param canonical_subsampling Oversampling factor per axis for high-res canonical grid. Default 8.
#' @param generate_plots Logical flag to trigger diagnostic multi-stage slice plots.
#' @param plot_filename Output filename for stage diagnostic plot.
#'
#' @return A list with true, hard, and mixel volumes, error percentages, and stage masks.
#' @export
sim_pve <- function(shape = "sphere", 
                    params = list(r = 4.0), 
                    res_vox = c(1.0, 1.0, 3.0), 
                    offset = c(0.0, 0.0, 0.0), 
                    rot_deg = c(0.0, 0.0, 0.0),
                    int_structure = 100, 
                    int_background = 50, 
                    noise_sd = 2.0,
                    canonical_subsampling = 8, 
                    generate_plots = FALSE, 
                    plot_filename = "mri_pve_stages.png") {
  
  # --- Verifications ---
  if (length(res_vox) != 3 || any(res_vox <= 0)) stop("res_vox must be a positive vector of length 3 (dx, dy, dz).")
  if (length(offset) != 3 || length(rot_deg) != 3) stop("offset and rot_deg must be length 3 vectors.")
  if (int_structure == int_background) stop("Structure and background intensities cannot be equal.")
  
  # =========================================================================
  # STEP 1: Generate Canonical Shape in Supersampled Isotropic Grid (0s & 1s)
  # =========================================================================
  buffer <- max(res_vox) * 4
  if (shape == "sphere") { max_dim <- params$r * 2
  } else if (shape == "ellipsoid") { max_dim <- max(params$a, params$b, params$c) * 2
  } else if (shape == "cylinder") { max_dim <- sqrt((2 * params$r)^2 + (params$h)^2) * 2
  } else if (shape == "ovoid") { max_dim <- max(params$a, params$b, params$c) * 2
  } else if (shape == "irregular") { max_dim <- (params$base_r + params$amp1 + params$amp2) * 2
  } else if (shape == "laminar") { max_dim <- max(c(params$thickness, params$curve_amp * 2)) * 2
  } else { stop("Unknown shape designated.") }

  total_span_mm <- max_dim + (2 * max(abs(offset))) + buffer
  
  # Fine isotropic resolution for high-fidelity analytical ground truth
  iso_res <- min(res_vox) / canonical_subsampling
  n_canon <- as.integer(ceiling(total_span_mm / iso_res))
  if (n_canon %% 2 == 0) n_canon <- n_canon + 1
  
  half_span <- (n_canon - 1) / 2
  canon_coords <- seq(-half_span, half_span) * iso_res
  
  # Grid coordinates centered at origin
  canon_grid <- expand.grid(x = canon_coords, y = canon_coords, z = canon_coords)
  
  # Transform coordinates: Translate and Rotate
  rot_rad <- rot_deg * pi / 180
  pts_shifted <- matrix(c(canon_grid$x - offset[1], 
                          canon_grid$y - offset[2], 
                          canon_grid$z - offset[3]), ncol = 3)
  pts_rotated <- rotate_coords(pts_shifted, -rot_rad[1], -rot_rad[2], -rot_rad[3])
  
  # Generate continuous binary mask (1 inside, 0 outside)
  inside_canon <- is_inside_shape(pts_rotated, shape, params)
  canonical_mask <- array(as.numeric(inside_canon), dim = c(n_canon, n_canon, n_canon))
  
  # High-fidelity Ground Truth Physical Volume
  voxel_vol_canon <- iso_res^3
  v_true <- sum(canonical_mask) * voxel_vol_canon
  if (v_true == 0) v_true <- 0.0001

  # =========================================================================
  # STEP 2: Resample Canonical Mask to Target Anisotropic Grid (Prior Mask)
  # =========================================================================
  n_vox <- as.integer(ceiling(total_span_mm / res_vox))
  n_vox <- ifelse(n_vox %% 2 == 0, n_vox + 1, n_vox)
  
  x_vox <- (seq_len(n_vox[1]) - (n_vox[1] + 1) / 2) * res_vox[1]
  y_vox <- (seq_len(n_vox[2]) - (n_vox[2] + 1) / 2) * res_vox[2]
  z_vox <- (seq_len(n_vox[3]) - (n_vox[3] + 1) / 2) * res_vox[3]
  
  # Downsample canonical mask via numerical integration (fractional occupancy)
  occupancy_map <- array(0.0, dim = n_vox)
  sub_x <- seq(-res_vox[1]/2, res_vox[1]/2, length.out = canonical_subsampling)
  sub_y <- seq(-res_vox[2]/2, res_vox[2]/2, length.out = canonical_subsampling)
  sub_z <- seq(-res_vox[3]/2, res_vox[3]/2, length.out = canonical_subsampling)
  sub_offset_grid <- as.matrix(expand.grid(dx = sub_x, dy = sub_y, dz = sub_z))
  n_sub <- nrow(sub_offset_grid)

  for (i in seq_len(n_vox[1])) {
    for (j in seq_len(n_vox[2])) {
      for (k in seq_len(n_vox[3])) {
        v_ctr <- c(x_vox[i], y_vox[j], z_vox[k])
        pts_v <- matrix(rep(v_ctr, each = n_sub), ncol = 3) + sub_offset_grid
        
        # Apply transformation
        pts_v_shift <- matrix(c(pts_v[,1] - offset[1], 
                                pts_v[,2] - offset[2], 
                                pts_v[,3] - offset[3]), ncol = 3)
        pts_v_rot <- rotate_coords(pts_v_shift, -rot_rad[1], -rot_rad[2], -rot_rad[3])
        occupancy_map[i, j, k] <- sum(is_inside_shape(pts_v_rot, shape, params)) / n_sub
      }
    }
  }
  
  # Initial binary anatomical prior mask (derived directly from resampled prior)
  initial_prior_mask <- occupancy_map >= 0.5
  voxel_vol_target <- prod(res_vox)

  # =========================================================================
  # STEP 3: Intensity & Noise Modeling
  # =========================================================================
  voxel_intensity_pure <- occupancy_map * int_structure + (1.0 - occupancy_map) * int_background
  noise_matrix <- array(rnorm(prod(n_vox), mean = 0, sd = noise_sd), dim = n_vox)
  voxel_intensity_observed <- voxel_intensity_pure + noise_matrix

  # =========================================================================
  # STEP 4: Adapt Boundary via Contiguous Intensity Expansion/Contraction
  # =========================================================================
  # Global threshold reference
  thresh <- (int_structure + int_background) / 2
  above_thresh <- if (int_structure > int_background) {
    voxel_intensity_observed >= thresh
  } else {
    voxel_intensity_observed <= thresh
  }

  # Establish 6-connected 3D neighborhood helper
  get_6_neighbors <- function(mask) {
    dim_m <- dim(mask)
    padded <- array(FALSE, dim = dim_m + 2)
    padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 2:(dim_m[3]+1)] <- mask
    
    dilated <- padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 2:(dim_m[3]+1)] |
               padded[1:(dim_m[1]),   2:(dim_m[2]+1), 2:(dim_m[3]+1)] |
               padded[3:(dim_m[1]+2), 2:(dim_m[2]+1), 2:(dim_m[3]+1)] |
               padded[2:(dim_m[1]+1), 1:(dim_m[2]),   2:(dim_m[3]+1)] |
               padded[2:(dim_m[1]+1), 3:(dim_m[2]+2), 2:(dim_m[3]+1)] |
               padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 1:(dim_m[3])]   |
               padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 3:(dim_m[3]+2)]
    
    eroded <-  padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 2:(dim_m[3]+1)] &
               padded[1:(dim_m[1]),   2:(dim_m[2]+1), 2:(dim_m[3]+1)] &
               padded[3:(dim_m[1]+2), 2:(dim_m[2]+1), 2:(dim_m[3]+1)] &
               padded[2:(dim_m[1]+1), 1:(dim_m[2]),   2:(dim_m[3]+1)] &
               padded[2:(dim_m[1]+1), 3:(dim_m[2]+2), 2:(dim_m[3]+1)] &
               padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 1:(dim_m[3])]   &
               padded[2:(dim_m[1]+1), 2:(dim_m[2]+1), 3:(dim_m[3]+2)]
    
    list(dilated = dilated, eroded = eroded)
  }

  adapted_mask <- initial_prior_mask

  # A. Contiguous Expansion: Grow into adjacent external voxels if above threshold
  neighbors <- get_6_neighbors(adapted_mask)
  external_boundary <- neighbors$dilated & (!adapted_mask)
  expansion_candidates <- external_boundary & above_thresh
  adapted_mask[expansion_candidates] <- TRUE

  # B. Contiguous Retraction: Pull back from adjacent internal edge voxels if sub-threshold
  neighbors_rev <- get_6_neighbors(adapted_mask)
  internal_boundary <- adapted_mask & (!neighbors_rev$eroded)
  retraction_candidates <- internal_boundary & (!above_thresh)
  adapted_mask[retraction_candidates] <- FALSE

  # =========================================================================
  # STEP 5: Hard & Mixel Volume Calculations (Contiguous Domain Masked)
  # =========================================================================
  # 1. Hard volume from locally adapted contiguous mask
  v_hard <- sum(adapted_mask) * voxel_vol_target
  err_hard <- ((v_hard - v_true) / v_true) * 100

  # 2. Raw linear mixel occupancy estimated from observed intensities
  raw_est_occupancy <- (voxel_intensity_observed - int_background) / (int_structure - int_background)
  raw_est_occupancy[raw_est_occupancy < 0] <- 0.0
  raw_est_occupancy[raw_est_occupancy > 1] <- 1.0

  # 3. Define the spatial evaluation domain for mixels:
  # Include the adapted contiguous structure mask + its 1-voxel 6-connected neighbor shell
  nb_adapted <- get_6_neighbors(adapted_mask)
  mixel_domain_mask <- nb_adapted$dilated  # Contiguous interior + contiguous boundary

  # 4. Mask the occupancy map so distant/non-contiguous background voxels are strictly zeroed
  mixel_occupancy_masked <- array(0.0, dim = n_vox)
  mixel_occupancy_masked[mixel_domain_mask] <- raw_est_occupancy[mixel_domain_mask]

  # 5. Calculate spatially restricted mixel volume
  v_mixel <- sum(mixel_occupancy_masked) * voxel_vol_target
  err_mixel <- ((v_mixel - v_true) / v_true) * 100

  # =========================================================================
  # STEP 6: Multi-Stage Diagnostic Plots (5 Rows x 3 Planes Layout)
  # =========================================================================
  if (generate_plots) {
    gray_pal <- gray.colors(256, start = 0, end = 1)
    heat_pal <- heat.colors(256, rev = TRUE) # Standard R heat palette (0=cool/yellow, 1=red/hot)
    # 1200 x 2000 Canvas: 5 Rows (Stages) x 3 Columns (Axial, Coronal, Sagittal)
    png(filename = plot_filename, width = 1200, height = 2000, res = 110)    
    # 5 rows, 3 columns layout matrix
    layout(matrix(1:15, nrow = 5, ncol = 3, byrow = TRUE))
    par(mar = c(3.5, 3.5, 2.5, 1.0), xaxs = "i", yaxs = "i")
    # Helper function to render a single 3-plane row for a given 3D matrix
    render_stage_row <- function(vol_data, x_coords, y_coords, z_coords, stage_title, 
                                 overlay_mask = NULL, palette = gray_pal, zlim = NULL) {
      mid_x <- round(length(x_coords) / 2)
      mid_y <- round(length(y_coords) / 2)
      mid_z <- round(length(z_coords) / 2)
      # 1. Axial plane (X vs Y)
      axial <- vol_data[, , mid_z]
      image(x_coords, y_coords, axial, col = palette, zlim = zlim,
            main = paste(stage_title, "- Axial"), xlab = "X (mm)", ylab = "Y (mm)", asp = 1)
      if (!is.null(overlay_mask)) {
        contour(x_coords, y_coords, overlay_mask[, , mid_z], levels = 0.5,
                col = "red", lwd = 2, add = TRUE, drawlabels = FALSE)
      }
      # 2. Coronal plane (X vs Z)
      coronal <- vol_data[, mid_y, ]
      image(x_coords, z_coords, coronal, col = palette, zlim = zlim,
            main = paste(stage_title, "- Coronal"), xlab = "X (mm)", ylab = "Z (mm)")
      if (!is.null(overlay_mask)) {
        contour(x_coords, z_coords, overlay_mask[, mid_y, ], levels = 0.5,
                col = "red", lwd = 2, add = TRUE, drawlabels = FALSE)
      }
      # 3. Sagittal plane (Y vs Z)
      sagittal <- vol_data[mid_x, , ]
      image(y_coords, z_coords, sagittal, col = palette, zlim = zlim,
            main = paste(stage_title, "- Sagittal"), xlab = "Y (mm)", ylab = "Z (mm)")
      if (!is.null(overlay_mask)) {
        contour(y_coords, z_coords, overlay_mask[mid_x, , ], levels = 0.5,
                col = "red", lwd = 2, add = TRUE, drawlabels = FALSE)
      }
    }

    # Row 1: Stage 1 - High-Res Canonical
    render_stage_row(canonical_mask, canon_coords, canon_coords, canon_coords, "1. Canonical High-Res")
    # Row 2: Stage 2 - Resampled Prior Mask
    render_stage_row(initial_prior_mask, x_vox, y_vox, z_vox, "2. Resampled Prior")
    # Row 3: Stage 3 - Intensity + Noise
    render_stage_row(voxel_intensity_observed, x_vox, y_vox, z_vox, "3. Intensity + Noise")
    # Row 4: Stage 4 - Adapted Boundary Overlay on Intensity
    render_stage_row(voxel_intensity_observed, x_vox, y_vox, z_vox, "4. Adapted Boundary", overlay_mask = adapted_mask)
    # Row 5: Stage 5 - Mixel Occupancy Heatmap (Restricted to Contiguous Domain)
    render_stage_row(mixel_occupancy_masked, x_vox, y_vox, z_vox, "5. Mixel Occupancy", 
                     palette = heat_pal, zlim = c(0, 1))

    dev.off()
  }

  list(
    true_volume = v_true,
    hard_volume = v_hard,
    hard_error_pct = err_hard,
    mixel_volume = v_mixel,
    mixel_error_pct = err_mixel,
    voxel_dims = n_vox,
    res_vox = res_vox
  )
}
