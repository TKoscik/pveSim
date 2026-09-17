\name{sim_pve}
\alias{sim_pve}
\title{Execute Advanced 3D MRI Partial Volume Error (PVE) Simulation}
\usage{
  sim_pve(shape = "sphere", params = list(r = 4), res_inplane = 1, res_z = 3,
          offset = c(0, 0, 0), rot_deg = c(0, 0, 0), int_structure = 100,
          int_background = 50, noise_sd = 2, subsamples = 4,
          generate_plots = FALSE, plot_filename = "mri_slice_output.png")
}
\arguments{
  \item{shape}{A character string defining the target structure. Must be one of
    \code{"sphere"}, \code{"ellipsoid"}, \code{"cylinder"}, \code{"egg"},
    \code{"irregular"}, or \code{"laminar"}.}
  
  \item{params}{A named list containing the physical size parameters matching the chosen shape.}
  
  \item{res_inplane}{A positive numeric value detailing the in-plane X and Y spatial resolution. Defaults to \code{1.0}.}
  
  \item{res_z}{A positive numeric value detailing the slice thickness (Z resolution). Defaults to \code{3.0}.}
  
  \item{offset}{A numeric vector of length 3 detailing the sub-voxel \code{c(x, y, z)} translation displacement offsets. Defaults to \code{c(0, 0, 0)}.}
  
  \item{rot_deg}{A numeric vector of length 3 detailing the roll, pitch, and yaw rotation angles in degrees. Defaults to \code{c(0, 0, 0)}.}
  
  \item{int_structure}{A numeric value defining the pure baseline MR signal intensity of the structure. Defaults to \code{100}.}
  
  \item{int_background}{A numeric value defining the pure baseline MR signal intensity of the background tissue. Defaults to \code{50}.}
  
  \item{noise_sd}{A non-negative numeric value defining the standard deviation of the ambient Gaussian electronic/thermal MR noise. Defaults to \code{2.0}.}
  
  \item{subsamples}{An integer value controlling the supersampling resolution (subdivisions per voxel axis) used for numerical fractional integration. Total sub-voxels evaluated per voxel is \code{subsamples^3}. Defaults to \code{4}.}
  
  \item{generate_plots}{A logical flag to manually trigger PNG slice generation through the spatial center of the voxel matrix array. Defaults to \code{FALSE}.}
  
  \item{plot_filename A character string specifying the file pathway where the diagnostic cross-section PNG layout should be exported if \code{generate_plots = TRUE}. Defaults to \code{"mri_slice_output.png"}.}
}
\value{
  A standard named list containing the following quantitative experimental metrics:
    \item{\code{true_volume}{The calculated baseline ground-truth continuous volume footprint.}
    \item{\code{hard_volume}{The estimated volume calculated via traditional binary 50\% threshold segmentation.}
    \item{\code{hard_error_pct}{The resulting percentage error under hard thresholding (positive values indicate overestimation, negative indicate underestimation).}
    \item{\code{mixel_volume}{The estimated volume calculated using soft fractional un-mixing (The Linear Mixel Model).}
    \item{\code{mixel_error_pct}{The resulting percentage error under soft fractional un-mixing.}
    \item{\code{boundary_voxels_count}{The absolute number of voxels holding a blended tissue profile.}
    \item{\code{contrast_to_noise_ratio}{The explicit Contrast-to-Noise Ratio (CNR) profile of the simulation run.}
    \item{\code{voxel_dimensions}{A vector tracking the dimensions of the generated simulation grid array.}
}
\description{
  Runs a comprehensive, multi-layered partial volume error simulation. The function maps
  continuous 3D geometries onto a discrete, customizable, anisotropic voxel grid, calculates
  mixed-tissue signal intensities based on structural and background contrast, layers on thermal
  MR noise, and returns both binary thresholded and linear mixel model volume measurements.
}
\examples{
  # Simulate an anisotropic acquisition of an ellipsoid structure with a 15-degree tilt
  my_sim <- sim_pve(
    shape = "ellipsoid", 
    params = list(a = 5.0, b = 3.5, c = 2.5),
    res_inplane = 1.0, res_z = 2.5,
    rot_deg = c(15, 0, 0),
    generate_plots = FALSE
  )
  print(my_sim$hard_error_pct)
}
