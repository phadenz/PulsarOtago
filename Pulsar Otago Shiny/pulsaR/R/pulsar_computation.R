#================================================================================
#     R implementation of PULSAR algorithm (Merriam & Wachter, 1982)

#     Copyright (C) 2020  Patricia Haden
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#     Merriam, G., & Wachter, K. (1982). Algorithms for the study of episodic
#     hormone secretion. The American journal of physiology, 243 4, E310-8.
#================================================================================

#///////////////////////////////////////////////////////////////////////////////
# Top level computation
#///////////////////////////////////////////////////////////////////////////////

#================================================================================
#' analyse_pulsar_data
#'
#' Top level Pulsar routine. Accepts prepared data frame, performs pulsar peak
#' location analysis on all data series, by passing to function pulsar_one. Returns
#' list of analysis output data structures.
#'
#' @param  hormone_df Pulsar format data frame as produced by read_pulsar_input_file
#' or read_pulsar_input_folder file. Column 1 = sample number; column 2 = time, columns
#' 3 to n are individual data series.
#'
#' @param smoothing_fraction Proportional width of the total time of session
#' to be used as the smoothing window for lowess smoothing. Should be between 0 and 1
#'
#' @param g_values Five-item vector of double defining the cutoff function g(n)
#' for peak identification (cf. Merriam and Wachter, 1982)
#'
#' @param extinction_threshold Minimum value allowed in data series. Non-negative
#' values less than extinction_threshold are replaced with extinction_threshold in
#' function clean_data
#'
#' @param peak_split_depth Size of drop and subsequent rise (in rescaled concentration
#' units) within a first pass peak that causes the peak to be redefined as two overlapping
#' events. See fsm_split_peaks for computational details.
#'
#' @param sdr_coef Three-item vector of double defining the coefficients in the equation
#' used to estimate assay variability for rescaling (see function gen_rescaled for details).
#' Terms should be in the order: quadratic, linear, constant.
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @param sdr_assay Constant estimate of assay variation. Not used in the current implementation
#'
#' @param n_steps Controls iteration of lowess algorithm. Default value of 6 replicates modified
#' lowess of Merriam and Wachter, 1982.
#'
#' @param y_axis_max If NULL (default) plot y-axis limits are set to the min/max of their
#' individual data series. If a value is provided, y-axis limits are set to (0, provided value).
#'
#' @return List of outputs from pulsar_one. Each containts fields plot, peak_descriptors,
#' peak_summary and animal_summary.
#'
#' @export
#================================================================================
analyse_pulsar_data <- function(hormone_df,
                                smoothing_fraction,
                                g_values,
                                extinction_threshold,
                                peak_split_depth,
                                sdr_coef,
                                nearest_nadir_distance,
                                sdr_assay = NULL,
                                n_steps = const$n_steps_constant,
                                y_axis_max = NULL)
{
  # Generate warnings if any suspicious data features are found (e.g. time values
  # not monotonically increasing)
  check_data_integrity(hormone_df)

  # Save the param values every time called. We'll load the last used values at
  # start up in the shiny app

  # generates param_list.Rdata
  save_parameters(smoothing_fraction,
                  g_values,
                  extinction_threshold,
                  peak_split_depth,
                  sdr_coef,
                  nearest_nadir_distance)

  # Parse df. Rule is that first two cols are count and time (see pulsar_constants.R)
  # The remaining columns are animal data, one animal per column
  n_animals <- ncol(hormone_df) - const$n_utility_columns
  sample_times <- as.numeric(hormone_df[,2])

  # Prepare loop driver sequence from 3 to n, the animal columns
  animal_column_indices <- (const$n_utility_columns + 1):(const$n_utility_columns + n_animals)

  # Prepare storage structure to gather individual animal output structures
  pulsar_outputs <- list()
  output_count <- 1

  # Process each animal
  for (a in animal_column_indices)
  {
    # Column names in the prepared input df should be the animal ids
    # (see read_input_file)

    # Pull the current id
    animal_id <- colnames(hormone_df)[a]

    # Pull the current concentration values into a vector
    animal_conc <- as.numeric(hormone_df[,a])

    # Combine time (common to all animals) and current concentration values
    # into a df
    animal_df <- data.frame(time = sample_times, conc = animal_conc)

    # Pass to the main computational routine pulsar_one. Optional argument
    # sdr_assay can be overriden for fised error term (as opposed to interpolated)
    animal_output <- pulsar_one(animal_id,
                                animal_df,
                                smoothing_fraction,
                                g_values,
                                extinction_threshold,
                                peak_split_depth,
                                sdr_coef,
                                nearest_nadir_distance,
                                y_axis_max = y_axis_max)


    # pulsar_one returns a list comprising one ggplot and two dfs, individual
    # peak data and peak summary. We build a list of these lists

    # NB: Use [[]] to build list of lists
    pulsar_outputs[[output_count]] <- animal_output
    output_count <- output_count + 1

  } # end for all animal columns

  # Return list of lists containing results for each animal
  return(pulsar_outputs)
}

#================================================================================
#' pulsar_one
#'
#' Accepts prepared data frame, with single data series. Performs pulsar peak
#' location analysis on single data series
#'
#' @param animal_id Unique string identifier for data series
#'
#' @param  animal_df Data frame. Column 1 = time; column 2 = data (hormone readings)
#'
#' @param smoothing_fraction Proportional width of the total time of session
#' to be used as the smoothing window for lowess smoothing. Should be between 0 and 1
#'
#' @param g_values Five-item vector of double defining the cutoff function g(n)
#' for peak identification (cf. Merriam and Wachter, 1982)
#'
#' @param extinction_threshold Minimum value allowed in data series. Non-negative
#' values less than extinction_threshold are replaced with extinction_threshold in
#' function clean_data
#'
#' @param peak_split_depth Size of drop and subsequent rise (in rescaled concentration
#' units) within a first pass peak that causes the peak to be redefined as two overlapping
#' events. See fsm_split_peaks for computational details.
#'
#' @param sdr_coef Three-item vector of double defining the coefficients in the equation
#' used to estimate assay variability for rescaling (see function gen_rescaled for details).
#' Terms should be in the order: quadratic, linear, constant.
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @param sdr_assay Constant estimate of assay variation. Not used in the current implementation
#'
#' @param n_steps Controls iteration of lowess algorithm. Default value of 6 replicates modified
#' lowess of Merriam and Wachter, 1982.
#'
#' @param y_axis_max If NULL (default) plot y-axis limits are set to the min/max of their
#' individual data series. If a value is provided, y-axis limits are set to (0, provided value).
#'
#' @return Named list containing fields plot, peak_descriptors,
#' peak_summary and animal_summary.
#'
#' @export
#================================================================================

pulsar_one <- function(animal_id,
                        animal_df,
                        smoothing_fraction,
                        g_values,
                        extinction_threshold,
                        peak_split_depth,
                        sdr_coef,
                        nearest_nadir_distance,
                        sdr_assay = NULL,
                        n_steps = const$n_steps_constant,
                        y_axis_max = NULL)
{
  # !!! Data format must be time in col 1 and conc in col 2. Function check_data_integrity will
  # find some problems, but no guarantee.
  # Set column names to prepare for cleaning and plotting
  colnames(animal_df) = c("Time", "Concentration")

  # Remove rows with negative values or NA. Replace concentrations below threshold
  animal_df <- clean_data(animal_df, extinction_threshold)

  # Pull the vectors for easier argument passing
  times <- animal_df$Time
  raw_concentrations <- animal_df$Concentration

  # Call lowess with the provided parameters to get the smoothed concentrations
  smoothed_concentrations <- gen_smoothed(times, raw_concentrations, smoothing_fraction, n_steps)

  # Rescaled concentrations are raw-smooth/assay noise. The noise term can be a provided
  # constant, or a quadratic interpolation using provided coefficients.

  if (is.null(sdr_assay)) {
    rescaled_conc <- gen_rescaled(raw_concentrations, smoothed_concentrations, sdr_coef = sdr_coef)
  } else {
    rescaled_conc <- gen_rescaled(raw_concentrations, smoothed_concentrations, sdr_assay = sdr_assay)
  }

  # Following Merriam and Wachter (1982), a peak of n points is a consecutive
  # sequence of rescaled concentration values >= g(n), with lower values of g (i.e. the
  # cutoff criterion) for larger n. From each point in the series
  # we check for sequences of decreasing length from maximum n down to 1.
  # Function gen_peaks returns a list holding 0/1 for out/in a peak for each
  # point in the time series in one vector, and the length of the found sequence
  # in a second vector
  peak_points <- find_peaks(rescaled_conc, g_values)

  # Peak sequences that dip dramatically then rise again are redefined as separate
  # overlapping peaks. Returns list of flags and scores for each point.
  # Start of subpeak marked in flag as -1.
  split_peak_points <- fsm_split_peaks(peak_points, peak_split_depth, rescaled_conc)

  # The preceding two calls assess concentration deltas and overlap to identify
  # points in peaks. When the first or last point is included in a peak it
  # may be ambiguous, as the concentration before (after) was not measured
  # making it impossible to determine the maximum point on the associated pulse
  # and thus the peak amplitude. In this case, we wish to exclude the point from
  # subsequent summary analysis, but still retain the knowledge that it was
  # highlighted as an event extreme from baseline. Function trim_end_peaks
  # determines whether first and last point peaks can be retained. The rule
  # is that, to be considered real peaks, they must either be followed (first
  # point) or preceded (last point) by a more extreme point, assuring that no
  # unmeasured higher value has occurred outside of the sample period.
  trimmed_peak_points <- trim_end_peaks(split_peak_points, rescaled_conc)

  # Function gen_peak_features takes the split peak list and walks it, returning a
  # data frame with a row for each identified peak containing first point, max
  # point and peak length
  peak_features <- gen_peak_features(trimmed_peak_points, rescaled_conc)

  # If you found any peaks, process them. Otherwise leave these elements NULL
  peak_descriptors <- NULL
  peak_summary <- NULL

  if (nrow(peak_features) > 0){
    #-------------------------------
    # Make individual peak descriptors  max point, max conc, amplitude.
    peak_descriptors <- gen_peak_descriptors(animal_id,
                                             times,
                                             raw_concentrations,
                                             smoothed_concentrations,
                                             peak_features,
                                             nearest_nadir_distance)


    #------------------------
    # Make peak summary table using individual peak features extracted above
    peak_summary <- gen_peak_summary(animal_id,
                                     times,
                                     raw_concentrations,
                                     smoothed_concentrations,
                                     peak_features,
                                     nearest_nadir_distance)
  }


  #------------------------
  # Make peak plot with peaks, ambiguous trim points, and nadirs coloured in
  peak_plot <- gen_peak_plot(animal_id,
                             animal_df,
                             smoothed_concentrations,
                             trimmed_peak_points,
                             peak_features,
                             nearest_nadir_distance,
                             y_axis_max)


  # Function summarise_animals computes summary statistics for the current animal's
  # data series
  animal_summary <- summarise_animal(animal_id, raw_concentrations)


  # Gather everything up for return
  pulsar_one_output <- list(peak_plot = peak_plot,
                            peak_descriptors = peak_descriptors,
                            peak_summary = peak_summary,
                            animal_summary = animal_summary)

  return(pulsar_one_output)

} # end pulsar_one




#================================================================================
#' summarise_animal
#'
#' Accepts vector of raw concentrations. Currently returns mean and stdev of raw
#' concentrations. Can be expanded in future as required.
#'
#' @param animal_id Unique string identifier for data series
#'
#' @param  raw_concentrations Vector of raw concentration values for single animal.
#'
#' @return Named list containing fields AnimalID, AvgConcentration and SdConcentration
#'
#' @import stats
#'
#' @export
#================================================================================
summarise_animal <- function(animal_id, raw_concentrations)
{
  avg_concentration <- mean(raw_concentrations)
  sd_concentration <- stats::sd(raw_concentrations)

  animal_summary <- data.frame(Animal_ID = animal_id,
                               Avg_Concentration = avg_concentration,
                               SD_Concentration = sd_concentration)
  return(animal_summary)
}


#================================================================================
#' gen_peak_plot
#'
#' Accepts intermediate results of a pulsar analysis of a single animal data
#' series (id, raw data df, smoothed baseline values, peak flags and peak
#' features) and produces three standardised output artefacts describing the
#' identified peaks/pulses. These are a ggplot with peak points marked, a named
#' list containing descriptors (first point, maximum concentration, amplitude,
#' etc.) for each peak, and a named list comtaining summary descritpors (e.g.
#' average amplitude, average peak interval, peak frequency) across all
#' identified peaks.
#'
#'
#' @param animal_id String serving as unique identified for the current data series
#'
#' @param animal_df A pulsar format data frame (Sample, Time, Concentration) containing the
#' raw data
#'
#' @param smoothed_concentrations The values of the corresponding baseline (lowess smoothed) as produced
#' by function gen_smoothed
#'
#' @param peak_points A named list containing peak_flags and peak_score defining the state of each
#' point in the data series, as produced by the peak locating functions find_peaks, fsm_split_peaks
#' and trim_end_peaks
#'
#' @param peak_features A named list of first and max points for each identified peak/pulse as returned
#' by function gen_peak_features
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @param y_axis_max If NULL (default) plot y-axis limits are set to the min/max of their
#' individual data series. If a value is provided, y-axis limits are set to (0, provided value).
#'
#' @return A named list of the three output artefacts plot, peak_descriptors and
#' peak_summary
#'
#' @import ggplot2
#' @export
#================================================================================
gen_peak_plot <- function(animal_id,
                          animal_df,
                          smoothed_concentrations,
                          peak_points,
                          peak_features,
                          nearest_nadir_distance,
                          y_axis_max = NULL)
{
  time_column_label <- "Time"
  concentration_column_label <- "Concentration"
  colnames(animal_df) = c(time_column_label, concentration_column_label)

  # Pull the vectors for tidier code
  times <- animal_df[[time_column_label]]
  raw_concentrations <- animal_df[[concentration_column_label]]


  # Generate the basic plot, before any peak information is added. This shows the
  # raw concentrations and smoothed (baseline) concentrations
  peak_plot <- gen_base_pulsar_plot(animal_df, animal_id,
                                    smoothed_concentrations,
                                    time_column_label, concentration_column_label,
                                    const$pulsar_font,
                                    y_axis_max)


  # If no peaks were found, just return the basic plot
  if (nrow(peak_features) == 0) {
    return(peak_plot)
  }


  # Add the peak points
  peak_flags <- peak_points$flag

  # All points that were included in a peak will be coloured red
  # Get x and y values for those points and store in a df
  peak_times <- times[peak_flags == const$flag_in_peak]
  peak_raw_concentrations <- raw_concentrations[peak_flags == const$flag_in_peak]
  peak_display_df = data.frame(peak_times, peak_raw_concentrations)



  # All points in ambiguous peaks will be coloured green
  ambig_peak_times <- times[peak_flags == const$flag_ambiguous_peak_point]
  ambig_peak_raw_concentrations <- raw_concentrations[peak_flags == const$flag_ambiguous_peak_point]
  ambig_peak_display_df = data.frame(ambig_peak_times, ambig_peak_raw_concentrations)

  # All near_nadir points will be coloured purple
  near_nadir_indices <- get_near_nadir_indices(raw_concentrations, peak_features, nearest_nadir_distance)
  near_nadir_times <- times[near_nadir_indices]
  near_nadir_raw_concentrations <- raw_concentrations[near_nadir_indices]
  near_nadir_display_df <- data.frame(near_nadir_times, near_nadir_raw_concentrations)


  peak_plot <- peak_plot +
    ggplot2::geom_point(data = peak_display_df,
                        ggplot2::aes(x = peak_times, y = peak_raw_concentrations),
                        color = "#d55e00", size = 4) +
    ggplot2::geom_point(data = ambig_peak_display_df,
                        ggplot2::aes(x = ambig_peak_times, y = ambig_peak_raw_concentrations),
                        color = "#f0e442", size = 4) +
    ggplot2::geom_point(data = near_nadir_display_df,
                        ggplot2::aes(x = near_nadir_times, y = near_nadir_raw_concentrations),
                        color = "#0072b2", size = 4)


  return(peak_plot)
}



#///////////////////////////////////////////////////////////////////////////////
# Input data preparatory processing
#///////////////////////////////////////////////////////////////////////////////

#================================================================================
#' gen_smoothed
#'
#' Calls lowess on time x concentration to generate smoothing line. The lowess method is in
#' library stats
#'
#' @param times Vector of time values
#'
#' @param  raw_concentrations Vector of raw concentration values for single animal.
#'
#' @param smoothing_fraction Proportional width of the total time of session
#' to be used as the smoothing window for lowess smoothing. Should be between 0 and 1
#'
#' @param n_steps Controls iteration of lowess algorithm. Value of 6 replicates modified
#' lowess of Merriam and Wachter, 1982.
#'
#' @return Vector of smoothed (i.e. baseline) values.
#'
#' @import stats
#' @export
#================================================================================
gen_smoothed <- function(times,
                         raw_concentrations,
                         smoothing_fraction,
                         n_steps = const$n_steps_constant)
{

  # Call lowess
  # lowess(x, y = NULL, f = 2/3, iter = 3, delta = 0.01 * diff(range(x)))

  smoothed <- stats::lowess(times, raw_concentrations, f = smoothing_fraction, iter = n_steps)

  # Return the smoothed y values
  smoothed_concentrations <- smoothed$y
  return(smoothed_concentrations)
}

#================================================================================
#' gen_rescaled
#'
#' Normalises raw concentrations as (raw-smooth)/noise.
#' Noise can be either a fixed estimate (if sdr_assay is supplied; not used in this
#' implementation) or an interpolated value computed for each data point as ax^2 + bx + c, if
#' sdr_coef(a,b,c) is supplied. See Merriam and Wachter, 1982 for how to compute
#' of a, b and c from assay replicates.
#'
#'
#' @param  raw_concentrations Vector of raw concentration values for single animal.
#'
#' @param  smoothed_concentrations Vector of computed baseline values. In current implementation
#' these are the output of function gen_smoothed, but alternative baseline computations could
#' be employed.
#'
#' @param sdr_coef Three-item vector of double defining the coefficients in the equation
#' used to estimate assay variability for rescaling (see function gen_rescaled for details).
#' Terms should be in the order: quadratic, linear, constant.
#'
#' @param sdr_assay Constant estimate of assay variation. Not used in the current implementation
#'
#' @return Vector of nomrmalise values.
#'
#' @import stats
#' @export
#================================================================================
gen_rescaled <- function(raw_concentrations,
                         smoothed_concentrations,
                         sdr_coef,
                         sdr_assay = NULL)
{

  if (length(raw_concentrations) != length(smoothed_concentrations))
  {
    stop("Different numbers of raw and baseline concentrations.")
  }

  n_observations <- length(raw_concentrations)

  # Use either a supplied constant or a quadratic estimate from the supplied
  # coefficients

  if (!is.null(sdr_assay)){
    sdr <- sdr_assay
    sdr_values <- rep(sdr, n_observations)
  } else {
    sdr_quad = sdr_coef[1]
    sdr_linear = sdr_coef[2]
    sdr_const = sdr_coef[3]

    sdr_values <- (sdr_quad * (raw_concentrations^2) +
                     sdr_linear * raw_concentrations +
                     sdr_const) * 0.01
  }

  # Either way, rescale as (raw - smooth)/sdr
  rescaled <- (raw_concentrations - smoothed_concentrations)/sdr_values

  return(rescaled)
}



#///////////////////////////////////////////////////////////////////////////////
# Peak Identification
#///////////////////////////////////////////////////////////////////////////////

#================================================================================
#' find_peaks
#'
#' Using the g_values for function g(n), a peak is any set of n
#' consecutive points whose rescaled values are greater than g(n). Function g(n) decreases
#' with increasing n, so you will identify narrower high peaks and wider low peaks.
#'
#'
#' @param  rescaled_conc Vector of rescaled concentration values as returned from gen_rescaled
#' or equivalent.
#'
#' @param  g_values Numeric values of function g(n). Values are rescaled cutoffs for sequence of length
#' n. In this implementation, n ranges from 1 to 5, inclusive.
#'
#' @return Named list containing peak_flags and peak_score.
#' \itemize{
#' \item peak_flags: For each point in input vector, 0 if point is not in a
#' peak/pulse as defined by g(n) and 1 if point is in peak/pulse.
#' \item peak_score: For each point in input vector contains the length (number
#' of consecutive points) in the peak containing the point, or 0 if the point
#' is not identified as being in a peak.
#' }
#'
#' @export
#================================================================================
find_peaks <- function(rescaled_conc, g_values)
{
  # Prepare
  n_observations <- length(rescaled_conc)
  peak_flag <- rep(0, n_observations)
  peak_score <- rep(0, n_observations)

  # For each point in the rescaled concentration vector, check to see
  # if it is the start of a sequence of qualifying points for lengths
  # 1 to 5 (the number of g-values). When a sequence is identified, mark
  # the corresponding elements of peak_flag and peak_score, and jump
  # ahead to the next point not contained in the sequence

  # Classic array-walking pattern
  i = 0

  # Walk the whole vector
  while (i <= n_observations)
  {
    i <- i+1

    # Check for sequence from 5 down to 1
    for (g_n in seq(5,1,-1))
    {
      cutoff <- g_values[g_n]

      points_found <- 0
      points_target <- g_n

      peak_walker <- i

      while ((peak_walker <= n_observations) &
             (rescaled_conc[peak_walker] >= cutoff))
      {
        points_found <- points_found + 1
        peak_walker <- peak_walker + 1
      }

      # Move back to the last good point
      peak_walker <- peak_walker - 1

      if (points_found >= points_target)
      {
        for (point_marker in i:peak_walker)
        {
          peak_flag[point_marker] <- 1
          peak_score[point_marker] <- points_found
        }
        i <- peak_walker
        break;
      } # end if peak found
    } #end for each g(n)
  } # end for each point

  # Bundle up flag and score (length) vectors for convenient return
  peak_points <- list(flag = peak_flag,
                      score = peak_score)
  return(peak_points)

}

#================================================================================
#' fsm_split_peaks
#'
#' Walks the identified peak points (from function find_peaks) looking for overlapping
#' peaks that have incorrectly been identified as members of a single sequence.
#' Following Merriam (1982) we decide this has occurred if, within a sequence,
#' the rescaled values drop by parameter peak_split_depth, and then rise by
#' peak_split_depth. The break is defined as the low point of that dip and return.
#' This implementation uses a four-state fsm rather than the linear algorithm of the original PULSAR
#' code.
#'
#' @param  peak_points Named list containing peak_flags and peak_scores as returned by function
#' find_peaks
#'
#' @param peak_split_depth Double defining the size (down and up) of change in the rescaled concentration
#' values of a peak sequence that triggers separation into separate peaks.
#'
#' @param  rescaled_conc Vector of normalised distances from baseline, as returned by function gen_rescaled
#' or equivalent.
#'
#' @return Named list containing peak_flags and peak_score.
#' \itemize{
#' \item peak_flags: For each point in input vector, 0 if point is not in a
#' peak/pulse as defined by g(n) and 1 if point is in peak/pulse.
#' \item peak_score: For each point in input vector contains the length (number
#' of consecutive points) in the peak containing the point, or 0 if the point
#' is not identified as being in a peak.
#' }
#'
#' @export
#================================================================================
fsm_split_peaks <- function(peak_points, peak_split_depth, rescaled_conc)
{

  # For readability
  st_seeking <- 0
  st_climbing <- 1
  st_split_check <- 2
  st_splitting <- 3

  main_walker <- 1
  current_peak_loc <-1
  current_peak_conc <- rescaled_conc[current_peak_loc]
  current_dip_loc <- 1
  current_dip_conc <- rescaled_conc[current_dip_loc]


  peak_flags <- peak_points$flag
  n_observations <- length(peak_flags)
  current_state <- st_seeking

  # Update flags
  while(main_walker <= n_observations)
  {
    next_loc <- main_walker + 1
    if (next_loc > n_observations) {
      # last flag stays unchanged
      break
    }

    current_conc <- rescaled_conc[main_walker]
    next_conc <- rescaled_conc[next_loc]
    delta_conc <- next_conc - current_conc
    next_flag <- peak_flags[next_loc]


    if (current_state == st_seeking) {
      if (next_flag == const$flag_in_peak){
        current_state <- st_climbing
      }

    } else if (current_state == st_climbing){
      if (delta_conc > 0) {
        if (next_conc > current_peak_conc){
          current_peak_loc <- next_loc
          current_peak_conc <- rescaled_conc[current_peak_loc]
        }
      } else if ((current_peak_conc - next_conc) > peak_split_depth) {
        current_state = st_split_check
        current_dip_loc <- next_loc
        current_dip_conc <- rescaled_conc[current_dip_loc]
      }
    } else if (current_state == st_split_check){
      if (delta_conc <= 0) {
        if (next_conc < current_dip_conc){
          current_dip_loc <- next_loc
          current_dip_conc <- rescaled_conc[current_dip_loc]
        }
      } else if ((current_conc - current_dip_conc) > peak_split_depth) {
        current_state <- st_splitting
      }
    } else if (current_state == st_splitting) {
      peak_flags[current_dip_loc] <- const$flag_overlap_peak_point # next pass, start new peak here
      main_walker <- current_dip_loc - 1 # start seeking from next loc after increment
      current_state = st_seeking
    }

    main_walker <- main_walker + 1
    if (peak_flags[main_walker] == 0) {
      current_state <- st_seeking
    }
  } # while main_walker < n

  # Bundle up flag and score (length) vectors for convenient return
  # !!! Score not currently used, but we're hanging on to it just in case
  peak_points <- list(flag = peak_flags,
                      score = peak_points$score)
  return (peak_points)
} # end fsm_split_peaks


#================================================================================
#' trim_end_peaks
#'
#' Ambiguous end points and their neighbour peaks are marked following Hackwell
#' rule: Is the next neighbour larger, and also still far enough from baseline
#' to qualify as a peak? If so, the maximum of this peak group is unambiguous,
#' and no change is needed. If not, every point in the peak needs to be flagged
#' as ambiguous. Ambiguous peak points are coloured differently and excluded
#' from the peak summary computation (see summarise_pulsar_result).
#'
#' @param  split_peak_points Named list containing peak_flags and peak_scores as returned by function
#' fsm_split_peaks
#'
#' @param  rescaled_conc Vector of normalised distances from baseline, as returned by function gen_rescaled
#' or equivalent.
#'
#' @return Named list containing peak_flags and peak_score.
#' \itemize{
#' \item peak_flags: For each point in input vector, 0 if point is not in a
#' peak/pulse as defined by g(n), 1 if point is in peak/pulse and 2 if the point
#' was originally in a peak, but the maximum of that peak cannot be unambiguously
#' identified
#' \item peak_score: For each point in input vector contains the length (number
#' of consecutive points) in the peak containing the point, or 0 if the point
#' is not identified as being in a peak.
#' }
#'
#' @export
#================================================================================
trim_end_peaks <- function(split_peak_points, rescaled_conc)
{
  peak_flags <- split_peak_points$flag
  n_observations <- length(peak_flags)

  # Check for 1st point flagged as peak. Logically cannot be overlap
  # restart, as those are only found in a dip
  if (peak_flags[1] == const$flag_in_peak){

      # Hackwell rule:
      # Is the next neighbour larger, and also still far enough from
      # baseline to qualify as a peak? If so, the maximum of this
      # peak group is unambiguous, and no change is needed. If not,
      # every point in the peak needs to be flagged as ambiguous. Ambiguous
      # peak points are coloured differently and excluded from the peak
      # summary computation (see other functions)

      first_conc <- rescaled_conc[1]
      next_neighbour_conc <- rescaled_conc[2]

      # Check for violations of the requirements
      if ((next_neighbour_conc < first_conc) | (peak_flags[2] != const$flag_in_peak)){

        # Mark this whole peak group as ambiguous (and make sure you don't
        # fall off the end in a very odd data series...)
        i <- 1
        while ((peak_flags[i] == const$flag_in_peak) & (i <= n_observations))
        {
          peak_flags[i] <- const$flag_ambiguous_peak_point
          i <- i + 1
        }

      } # if first point is ambiguous

    } # if first point is flagged


    # The logic for the last point is the same, but the neighbour to check is
    # the previous point, and the group to be marked are before the ambiguous
    # last point.

    if (peak_flags[n_observations] == const$flag_in_peak){

      last_conc <- rescaled_conc[n_observations]
      prev_neighbour_conc <- rescaled_conc[n_observations - 1]

      # Check for violations. Last point is ambiguous is the previous
      # is not a peak, or if the previous is smaller than the last, even
      # if it is a peak
      if ((prev_neighbour_conc < last_conc) | (peak_flags[n_observations - 1] != const$flag_in_peak)) {

        # Mark this whole peak group as ambiguous (and make sure you don't
        # fall off the end in a very odd data series...)
        i <- n_observations
        while ((peak_flags[i] == const$flag_in_peak) & (i >= 1))
        {
          peak_flags[i] <- const$flag_ambiguous_peak_point
          i <- i - 1
        }
      } # last point is ambiguous
    } # if last point is flagged


    # Bundle up flag and score (length) vectors for convenient return
    # !!! Score not currently used, but we're hanging on to it just in case
    trimmed_peak_points <- list(flag = peak_flags,
                                score = split_peak_points$score)
    return (trimmed_peak_points)

} # end trim_ends


#================================================================================
#' gen_peak_features
#'
#' Accepts a peak_points list ($flag and $score) and produces
#' a descriptor df (first, max and length). The flags should be generated
#' by running find_peaks and fsm_split_peaks. Points not in a peak are 0.
#' Other peak flags show in peak, overlap point or ambiguous point as defined
#' by the corresponding constants (see file puslar_constants.R)
#'
#' @param  peak_points Named list containing peak_flags and peak_scores as returned by function
#' fsm_split_peaks and trim_end_peaks (if desired)
#'
#' @param  rescaled_conc Vector of normalised distances from baseline, as returned by function gen_rescaled
#' or equivalent.
#'
#' @return Data frame containing first point (FirstPoint), point of maximum
#'   concentration (MaxPoint) and length (Length) for each identified peak/pulse
#'   sequence. If no peaks are found, returns a df with 0 rows and 0 cols !!! TBD change to null?
#'
#' @export
#================================================================================

gen_peak_features <- function(peak_points, rescaled_conc)
{

  peak_flags <- peak_points$flag

  n_observations <- length(peak_flags)

  # R vectors are 1-indexed
  i <- 1
  n_peaks <- 1

  peak_first_point <- c()
  peak_max_point <- c()
  peak_length <- c()

  while (i <= n_observations)
  {
    # Find a peak
    while (((peak_flags[i] == 0) | (peak_flags[i] == const$flag_ambiguous_peak_point)) & (i <= n_observations))
    {
      i <- i + 1
    }

    if (i > n_observations) {break}

    peak_first_point[n_peaks] <- i

    # Count the points in this peak
    count <- 0
    while ((i <= n_observations) & (peak_flags[i] == const$flag_in_peak))
    {
      count <- count + 1
      i <- i + 1
    }

    # Store count
    peak_length[n_peaks] <- count

    # Get and store loc of max concentration
    first_loc <- peak_first_point[n_peaks]
    max_loc <- first_loc
    max_conc <- rescaled_conc[first_loc]

    for (p in first_loc:(first_loc + count - 1))
    {
      if (rescaled_conc[p] > rescaled_conc[max_loc])
      {
        max_loc <- p
      }
    }
    peak_max_point[n_peaks] <- max_loc


    # Prepare for next peak
    n_peaks <- n_peaks + 1

    if ((i <= n_observations) & (peak_flags[i] == const$flag_overlap_peak_point)) {
      peak_flags[i] <- 1
    }
  } # i <= n

  peak_df <- data.frame(FirstPoint = peak_first_point,
                        MaxPoint = peak_max_point,
                        Length = peak_length)
  return(peak_df)

} # end gen_peak_features


#===============================================================================
#' get_near_nadir_indices
#'
#' For each located peak/pulse, finds the index of the "nearest nadir" point. This is the
#' point with the lowest raw concentration in a window p samples prior to the maximum
#' concentration point of the peak. The window size p is a parameter of the algorithm passed
#' as argument nearest_nadir_distance. These values are used in the computation of amplitude.
#'
#' @param raw_concentrations Raw concentration values
#'
#' @param peak_features A named list of first and max points for each identified peak/pulse as returned
#' by function gen_peak_features
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @return Vector of indices of the nearest nadir of each peak.
#'
#' @export
#===============================================================================
get_near_nadir_indices <- function(raw_concentrations, peak_features, nearest_nadir_distance)
{

  near_nadir_indices <- c()

  # This definition replaced with constant 05-02-21
  # find mean sequence length. These values are held in peak_features$Length
  # This is the distance to go back when looking for the nearest nadir.
  #mean_sequence_length <- mean(peak_features$Length)
  #mean_sequence_length <- ceiling(mean_sequence_length) # round up

  for (max_point_index in peak_features$MaxPoint)
  {
    # Find where to start looking for the nadir
    first_index <- max_point_index - nearest_nadir_distance
    if (first_index < 1) {
      first_index <- 1
    }

    # In the range from first index to max index, which value has the lowest raw concentration?
    near_nadir_index <- which.min(raw_concentrations[first_index:(max_point_index - 1)])

    # Translate that index value from the segment location to the whole data set location
    near_nadir_index <- near_nadir_index + first_index - 1

    near_nadir_indices <- c(near_nadir_indices, near_nadir_index)
  }

  return(near_nadir_indices)
}

#===============================================================================
#' get_near_nadir_conc
#'
#' Given the index of a peak maximum, returns the raw concentration at the "nearest nadir" point.
#' This is the point with the lowest raw concentration in a window p samples prior to the maximum
#' concentration point of the peak. The window size p is a parameter of the algorithm passed
#' as argument nearest_nadir_distance. These values are used in the computation of amplitude.
#'
#' @param raw_concentrations Raw concentration values
#'
#' @param max_point_index Integer value. Index of a point in the data vector for which the
#' nearest nadir concentration is to be calculated.
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @return Vector of indices of the nearest nadir of each peak.
#'
#' @export
#===============================================================================
get_near_nadir_conc <- function(raw_concentrations, max_point_index, nearest_nadir_distance)
{
  first_index <- max_point_index - nearest_nadir_distance
  if (first_index < 1) {
    first_index <- 1
  }
  near_nadir_conc <- min(raw_concentrations[first_index:(max_point_index - 1)])
  return(near_nadir_conc)
}

#================================================================================
#' gen_amplitude
#'
#' Generates amplitude for each found peak. Amplitude is computed as the
#' difference between the pulse's maximum raw concentration and the value of
#' the "nearest nadir". The nearest nadir is the lowest concentration value in
#' the previous nearest_nadir_distance points.
#' Amplitude is max raw conc minus raw conc at nearest nadir neighbour
#'
#' @param raw_concentrations Raw concentration values
#'
#' @param peak_features A named list of first and max points for each identified peak/pulse as returned
#' by function gen_peak_features
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the first point contained in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @return Vector of amplitudes for each peak
#'
#' @export
#================================================================================
gen_amplitude <- function(raw_concentrations, peak_features, nearest_nadir_distance)
{
  #return(raw_concentrations[peak_features$MaxPoint] - smoothed_concentrations[peak_features$MaxPoint])

  near_nadir_amplitudes <- c()

  # This definition replace with constant. 05-02-21
  # find mean sequence length. These values are held in peak_features$Length
  #mean_sequence_length <- mean(peak_features$Length)
  #mean_sequence_length <- ceiling(mean_sequence_length) # round up

  for (max_point_index in peak_features$MaxPoint)
  {
    # get concentration at maximum point
    max_conc <- raw_concentrations[max_point_index]

    # find lowest neighbour within nearest_nadir_distance previous points
    near_nadir_conc <- get_near_nadir_conc(raw_concentrations, max_point_index, nearest_nadir_distance)

    curr_amplitude <- max_conc - near_nadir_conc

    near_nadir_amplitudes <- c(near_nadir_amplitudes, curr_amplitude)
  }

  return(near_nadir_amplitudes)
}
#================================================================================
#' gen_peak_descriptors
#'
#' Generates summary information for found peaks using raw concentrations, smoothed concentrations
#' and peak_features (first and maximum points) as produced by method gen_peak_features. Summary
#' values are start times, maximum concentration time, maximum concentration
#' value, peak amplitude, and peak length for each identified peak/pulse.
#'
#' @param animal_id Unique string identifier for the data series
#'
#' @param times Time each assay sample was recorded
#'
#' @param raw_concentrations Raw concentration values
#'
#' @param smoothed_concentrations The values of the corresponding baseline, Currently lowess smoothed as produced
#' by function gen_smoothed
#'
#' @param peak_features A named list of first and max points for each identified peak/pulse as returned
#' by function gen_peak_features
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the first point contained in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @return Data frame containing peak descriptor values (columns) for each peak (rows)
#'
#' @export
#================================================================================
gen_peak_descriptors <- function(animal_id,
                                 times,
                                 raw_concentrations,
                                 smoothed_concentrations,
                                 peak_features,
                                 nearest_nadir_distance)
{
  # If there are no peak points returns a df with 0 rows
  if ((is.null(peak_features)) | (nrow(peak_features) == 0)){
    return(NULL)
  }

  # First column is animal_id for simplifying later analyses
  animal_id_vector <- rep(animal_id, nrow(peak_features))

  # Grab peak time and concentration features using info from final_peaks
  # final_peaks contains first point, max point and length for each id'd peak
  peak_start_time_vector <- times[peak_features$FirstPoint]
  peak_max_time_vector <- times[peak_features$MaxPoint]
  peak_max_conc_vector <- raw_concentrations[peak_features$MaxPoint]

  peak_amplitude_vector <- gen_amplitude(raw_concentrations, peak_features, nearest_nadir_distance)


  # Bundle individual peak features into handy data frame
  # peak_descriptor <- data.frame(AnimalID = animal_id,
  #                       StartTime = peak_start_time_vector,
  #                       MaxTime = peak_max_time_vector,
  #                       MaxConc = peak_max_conc_vector,
  #                       Amplitude = peak_amplitude_vector,
  #                       PeakLength = peak_features$Length)

  # 19-06-21 StarTime and PeakLength removed, as client doesn't want them
  peak_descriptor <- data.frame(Animal_ID = animal_id,
                                Peak_Time = peak_max_time_vector,
                                Max_Conc = peak_max_conc_vector,
                                Amplitude = peak_amplitude_vector)

  return(peak_descriptor)
}

#================================================================================
#' gen_peak_summary
#'
#' Generates summary information across all peaks: average amplitude and standard deviation
#' of amplitudes, average peak interval and standard deviation of peak intervals, and peak
#' frequency (peaks per hour)
#'
#' @param animal_id Unique string identifier for the data series
#'
#' @param times Time each assay sample was recorded
#'
#' @param raw_concentrations Raw concentration values
#'
#' @param smoothed_concentrations The values of the corresponding baseline, Currently lowess smoothed as produced
#' by function gen_smoothed
#'
#' @param peak_features A named list of first and max points for each identified peak/pulse as returned
#' by function gen_peak_features
#'
#' @param nearest_nadir_distance Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the first point contained in the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @return Data frame containing summary values across all peaks.
#'
#' @export
#================================================================================
gen_peak_summary <- function(animal_id,
                             times,
                             raw_concentrations,
                             smoothed_concentrations,
                             peak_features,
                             nearest_nadir_distance)
{

  # !!! TBD decide what happen to peak_features if there are no peak points
  # Currently returns a df with 0 rows, but maybe should return null?
  if ((is.null(peak_features)) | (nrow(peak_features) == 0)){
    return(NULL)
  }

  peak_amplitude_vector <- gen_amplitude(raw_concentrations, peak_features, nearest_nadir_distance)

  # Average and standard deviation of amplitude
  avg_amp <- mean(peak_amplitude_vector)
  sd_amp <- stats::sd(peak_amplitude_vector)


  # Identify durations of n_peaks-1 intervals (in t) between n_peaks peaks
  peak_intervals <- c()
  n_peaks <- nrow(peak_features)
  for (p in 1:(n_peaks-1))
  {
    firstPeakIndex = peak_features$MaxPoint[p]
    secondPeakIndex = peak_features$MaxPoint[p+1]
    interval <- times[secondPeakIndex] - times[firstPeakIndex]
    peak_intervals <- c(peak_intervals, interval)
  }

  # Summarise peak intervals
  avg_peak_interval <- mean(peak_intervals)
  sd_peak_interval <- stats::sd(peak_intervals)

  # Compute peak frequency as #peaks / (total_minutes/60)
  # for peaks/hour
  # !!! TBD CONFIRM. We define time as starting at the first reading here, not
  # at the first peak.

  total_minutes <- times[length(times)] - times[1]
  total_hours <- total_minutes/60
  peak_frequency <- n_peaks/total_hours

  # Bundle global peak summary values into handy data frame
  peak_summary <- data.frame(Animal_ID = animal_id,
                                Avg_Amp = avg_amp,
                                SD_Amp = sd_amp,
                                Avg_Peak_Interval = avg_peak_interval,
                                SD_Peak_Interval = sd_peak_interval,
                                Peak_Frequency_Hour = peak_frequency)
  return(peak_summary)

}
