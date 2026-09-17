# pveSim: MRI Partial Volume Error Simulator

`pveSim` is a lightweight, zero-dependency, base-R simulation package designed to quantify, model, and visualize the impact of **Partial Volume Effects (PVE)** and spatial alignment anomalies in 3D digital imaging, specifically tailored for Magnetic Resonance Imaging (MRI) contexts.

When continuous anatomical objects are sampled onto discrete voxel grids—especially anisotropic grids where slice thickness (Z) is larger than in-plane resolution (X, Y)—significant structural volume measurement bias is introduced. This package allows you to systematically model these geometric challenges, providing a rigorous sandbox to generate dummy data for downstream statistical pipelines, morphometry validation studies, and algorithmic threshold benchmarks.

## Key Architectural Strengths

* **Zero Tidyverse Dependencies:** Built purely using high-performance base-R linear algebra, arrays, matrices, and core numerical integration mechanics.
* **Complex Geometric Profiles:** Simulates traditional shapes (spheres, cylinders, ellipsoids), asymmetric deep nuclear bodies (ovoids/eggs), thin layered structures (cortical laminar ribbons), and highly irregular multi-frequency bumpy masses.
* **Realistic Scanner Realities:** Includes sub-voxel translation offsets (phase-shifting), multi-axial oblique rotations, independent axis anisotropy controls, variable tissue-to-background contrast settings, and Gaussian thermal/instrumentation noise parameters.
* **Dual Segmentation Outputs:** Contrasts traditional binary threshold classifications (`hard_error`) against advanced linear fractional mixture calculations (`mixel_error`).
* **Manual Diagnostic Export:** Generates rapid, native 3-plane (Axial, Coronal, Sagittal) cross-sectional PNG maps through the spatial core of your customized simulation space.

---

## Installation

You can install the development version of `pveSim` directly from your GitHub repository using `remotes` or `devtools`:

```R
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install pveSim from GitHub
remotes::install_github("TKoscik/pveSim")
```

---

## Core Functions Quick Reference

* `make_voxel_grid()`: Spits out a zero-centered 3D coordinate array mapping individual voxel center indices across custom anisotropic dimensions.
* `rotate_coords()`: Multiplies coordinate coordinates against a unified Roll-Pitch-Yaw mathematical matrix mapping rotation adjustments.
* `is_inside_shape()`: Vectorized intersection engine tracking whether 3D space vectors fall within or outside analytical shape thresholds.
* `sim_pve()`: Main wrapper pipeline controlling voxelization loops, noise application, segmentation estimation calculations, and PNG generation features.

---

## Usage Examples

### 1. Execute a Single Anisotropic Simulation Run
This code runs an anisotropic scan over an asymmetrical egg-shaped structure (simulating an irregular deep brain nucleus) tilted obliquely, and checks how hard thresholding handles the geometry:

```R
library(pveSim)

# Define an asymmetric egg shape parameter profile
egg_config <- list(a = 4.5, b = 3.5, c = 4.0, k = 0.15)

# Run simulation
result <- sim_pve(
  shape         = "egg",
  params        = egg_config,
  res_inplane   = 1.0,           # High in-plane resolution
  res_z         = 3.0,           # Thick clinical slices (3x Anisotropy)
  offset        = c(0.2, 0, 0),  # Sub-voxel shift
  rot_deg       = c(15, 12, 5),  # Random rotational tilt
  int_structure = 120,          # Clear white/gray matter intensity gaps
  int_background= 70,
  noise_sd      = 3.0            # Realistic instrument noise
)

# Examine directional error deviations
print(paste("Binary Thresholding Error Bias:", round(result$hard_error_pct, 2), "%"))
print(paste("Soft Mixel Model Error Bias:", round(result$mixel_error_pct, 2), "%"))
```

### 2. Generate a Manual Diagnostic Image Panel
By explicitly enabling `generate_plots = TRUE`, you can instantly output a 3-panel cross-section layout file detailing exactly how your structure interfaces with the noisy voxel columns:

```R
# Evaluate an irregular multi-frequency bumpy structure and output a slice file
irregular_config <- list(base_r = 4.0, amp1 = 0.8, freq1 = 4, amp2 = 0.4, freq2 = 2)

sim_pve(
  shape          = "irregular",
  params         = irregular_config,
  res_z          = 2.5,
  generate_plots = TRUE,
  plot_filename  = "irregular_structure_slice.png"
)
```

### 3. Loop Pipeline for Downstream Statistical Modeling
This script runs a quick Monte Carlo style batch experiment, varying rotation and anisotropy properties over a thin cortical laminar sheet structure, and applies a linear regression model to map out variance:

```R
set.seed(42)
iterations <- 30
laminar_config <- list(base_z = 0, curve_amp = 1.5, curve_freq = 0.5, thickness = 1.2)
dataset <- vector("list", iterations)

for(i in 1:iterations) {
  # Vary slice profiles (1x to 4x anisotropy layers) and angular head adjustments
  slice_thick <- sample(c(1.0, 2.0, 3.0, 4.0), 1)
  rand_rot    <- runif(3, min = 0, max = 45)
  rand_shift  <- runif(3, min = -0.5, max = 0.5)
  
  sim_run <- sim_pve(
    shape       = "laminar",
    params      = laminar_config,
    res_inplane = 1.0,
    res_z       = slice_thick,
    offset      = rand_shift,
    rot_deg     = rand_rot,
    noise_sd    = 2.5
  )
  
  dataset[[i]] <- data.frame(
    slice_thickness = slice_thick,
    max_tilt        = max(rand_rot),
    hard_error      = sim_run\$hard_error_pct,
    mixel_error     = sim_run\$mixel_error_pct
  )
}

# Aggregate list arrays into a standard base-R data frame
sim_df <- do.call(rbind, dataset)

# Fit an ordinary linear regression model to examine directional over/underestimation forces
pve_regression <- lm(abs(hard_error) ~ slice_thickness + max_tilt, data = sim_df)
summary(pve_regression)
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for complete details.
