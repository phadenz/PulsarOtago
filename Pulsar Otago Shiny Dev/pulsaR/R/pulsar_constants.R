const <- list()

# Flag states in peak_flags. Functions find_peaks,
# fsm_split_peak, trim_end_peak, get_peak_features and get_peak_descriptors
const$flag_overlap_peak_point = -1
const$flag_ambiguous_peak_point = -2
const$flag_in_peak = 1

# Appropriate value to match numeric results of modified
# lowess in Merriam and Wachter 1982
const$n_steps_constant = 6

# To support cross-platform good behaviour by ggplot.
const$pulsar_font <- "sans"


# Image size properties for combined image
const$mult_image_height = 1024
const$mult_image_width = 1280


# Columns 1 and 2 are sample number (just 1..n) and time
const$n_utility_columns = 2


# Default values based on Otago testing
# 19-06-21 Change to 0 for manuscript release
const$smoothing_fraction <- 1
const$g_values <- c(0,0,0,0,0)
const$extinction_threshold <- 0
const$peak_split_depth <- 0
const$sdr_coef = c(0,0,1)


const$sdr_assay = NULL # Not used in this implementation
const$n_steps = 6 # improves lowess accuracy over original value of 3

# Number of points to look back from a peak maximum to find the low point against
# which to compute amplitude. Typical LH data works well with 3
const$nearest_nadir_distance = 3
