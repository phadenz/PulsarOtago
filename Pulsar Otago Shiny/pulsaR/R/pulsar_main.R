# Top-level user-facing Pulsar methods run_pulsar and run_pulsar_batch
# for single command processing

#===================================================================================
#' run_pulsar
#'
#' Top-level user-facing pulsaR method.
#'
#' Accepts a file in Pulsar Otago input format (see below). Runs
#' complete analysis. Returns list of pulsar results (plots, peak descriptors,
#' peak summaries and animal summaries). Generates outfiles if
#' make_outfiles = TRUE (default).
#'
#'
#' Input file must be .csv.
#'
#' File format is:
#'
#' \itemize{
#' \item Line 1: Experiment name (any string) as the first element
#' \item Line 2: Column header line. Column headers can be any strings. Data
#' series column headers are treated as animal IDs, so should be unique.
#' \item Lines 3 to 500: Data values.
#'}
#' Columns should be:
#'
#' \itemize{
#' \item Column 1: Sample numbers (1 to n).
#' \item Column 2: Times (in minutes).
#' \item Columns 3 to n: Data values (typically hormone readings). Each data series (typically
#' animal) is in a single column.
#'}
#' All data series in a single file use the same time values.
#'
#' @param input_file_name Complete file path to a correctly formatted input file (see details below)
#'
#' @param param_list Named list of algorithm parameters as returned by function restore_parameters.
#' List can be built manually if desired, but field names must be as specified.
#'
#' \itemize{
#' \item smoothing_fraction: Proportional width of the total time of session
#' to be used as the smoothing window for lowess smoothing. Should be between 0 and 1.
#' \item g_values: Five-item vector of double defining the cutoff function g(n)
#' for peak identification (cf. Merriam and Wachter, 1982).
#' \item extinction_threshold: Minimum value allowed in data series. Non-negative
#' values less than extinction_threshold are replaced with extinction_threshold in
#' function clean_data.
#' \item peak_split_depth: Size of drop and subsequent rise (in rescaled concentration
#' units) within a first pass peak that causes the peak to be redefined as two overlapping
#' events. See fsm_split_peaks for computational details.
#' \item sdr_coef: Three-item vector of double defining the coefficients in the equation
#' used to estimate assay variability for rescaling (see function gen_rescaled for details).
#' Terms should be in the order: quadratic, linear, constant.
#' \item nearest_nadir_distance: Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point of the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#' \item sdr_assay: Constant assay variability estimate. Not used in this implementation; may
#' be null in the param_list.
#' \item n_steps: Lowess iteration term. Defaults to 6 if null in param_list. This value matches
#' numeric output of modified lowess as described by Merriam and Wachter, 1982.
#' }
#'
#' @param make_outfiles If TRUE (the default) output files are written to the working directory.
#'
#' @param zip_files If TRUE, output files are bundled into a single zip file.
#'
#' @param y_axis_max If NULL (default) plot y-axis limits are set to the min/max of their
#' individual data series. If a value is provided, y-axis limits are set to (0, provided value).
#'
#' @return List of pulsar result objects as produced by method analyse_pulsar_data. Each contains
#' plot, peak descriptors, peak summary and animal summary for one data series.
#'
#' @export
#===================================================================================
run_pulsar <- function(input_file_name, param_list = NULL, make_outfiles = TRUE, zip_files = FALSE, y_axis_max = NULL)
{
  # User can set up a param_list if they wish. Otherwise use last stored value,
  # or defaults, if no stored Rdata object
  if (is.null(param_list)) {
    param_list <- restore_param_list() # Last used or default
  }

  smoothing_fraction <- param_list[["smoothing_fraction"]]
  g_values <- param_list[["g_values"]]
  extinction_threshold <- param_list[["extinction_threshold"]]
  peak_split_depth <- param_list[["peak_split_depth"]]
  sdr_coef = param_list[["sdr_coef"]]
  nearest_nadir_distance = param_list[["nearest_nadir_distance"]]

  n_steps = param_list[["n_steps"]]
  if (is.null(n_steps)){
    n_steps = 6
    }
  sdr_assay <- param_list[["sdr_assay"]]

  # Returns list with ExperimentName and HormoneDF
  input_data <- read_pulsar_input_file(input_file_name)

  # If something goes wrong, read_pulsar_input_file will return null
  if (is.null(input_data)){
    warning(paste("Error processing data file", input_file_name))
    return(NULL)
  }

  # Pull components for cleaner code
  experiment_name <- input_data$ExperimentName
  hormone_df <- input_data$HormoneDF

  # Run analysis. Returns list of pulsar results: plot,peak descriptor list,
  # peak summary list, animal summary list
  all_results <- analyse_pulsar_data(hormone_df,
                                     smoothing_fraction,
                                     g_values,
                                     extinction_threshold,
                                     peak_split_depth,
                                     sdr_coef,
                                     nearest_nadir_distance,
                                     sdr_assay,
                                     n_steps,
                                     y_axis_max)


  # Write lists to usefully names files and generate png for each plot
  if (make_outfiles) {
    all_files <- write_multi_pulsar(experiment_name, param_list, all_results, zip_files = zip_files)
  }


  return (all_results)
}

#===================================================================================
#' run_pulsar_batch
#'
#' Top-level user-facing pulsaR method.
#'
#' Accepts a folder containing files in Pulsar Otago input format (see below).
#' All files must contain exactly one data series. Runs complete analysis.
#' Returns list of pulsar results. Generates outfiles if make_outfiles = TRUE
#' (default)
#'
#'
#' Input file must be .csv.
#'
#' File format is:
#'
#' \itemize{
#' \item Line 1: Experiment name (any string) as the first element.
#' \item Line 2: Column header line. Column headers can be any strings. Data
#' series column headers are treated as animal IDs, so should be unique.
#' \item Lines 3 to 500: Data values.
#'}
#' Columns should be:
#'
#' \itemize{
#' \item Column 1: Sample numbers (1 to n).
#' \item Column 2: Times (in minutes).
#' \item Columns 3 to n: Data values (typically hormone readings). Each data series (typically
#' animal) is in a single column.
#'}
#' All data series in a single file use the same time values.
#'
#' @param input_folder_name Complete path to a folder containing correctly formatted input files (see below)
#' each with exactly one data series.
#'
#' @param param_list Named list of algorithm parameters as returned by function restore_parameters.
#' List can be built manually if desired, but field names must be as specified.
#'
#' \itemize{
#' \item smoothing_fraction: Proportional width of the total time of session
#' to be used as the smoothing window for lowess smoothing. Should be between 0 and 1.
#' \item g_values: Five-item vector of double defining the cutoff function g(n)
#' for peak identification (cf. Merriam and Wachter, 1982).
#' \item extinction_threshold: Minimum value allowed in data series. Non-negative
#' values less than extinction_threshold are replaced with extinction_threshold in
#' function clean_data.
#' \item peak_split_depth: Size of drop and subsequent rise (in rescaled concentration
#' units) within a first pass peak that causes the peak to be redefined as two overlapping
#' events. See fsm_split_peaks for computational details.
#' \item sdr_coef: Three-item vector of double defining the coefficients in the equation
#' used to estimate assay variability for rescaling (see function gen_rescaled for details).
#' Terms should be in the order: quadratic, linear, constant.
#' \item nearest_nadir_distance: Amplitude is computed as the difference between the
#' maximum raw concentration in peak/pulse sequence and the raw concentration of the
#' "nearest nadir". The nearest nadir is the minimum raw concentration among n sample
#' points before the max point of the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#' \item sdr_assay: Constant assay variability estimate. Not used in this implementation; may
#' be null in the param_list.
#' \item n_steps: Lowess iteration term. Defaults to 6 if null in param_list. This value matches
#' numeric output of modified lowess as described by Merriam and Wachter, 1982.
#' }
#'
#' @param make_outfiles If TRUE (the default) output files are written to the working directory.
#'
#' @param y_axis_max If NULL (default) plot y-axis limits are set to the min/max of their
#' individual data series. If a value is provided, y-axis limits are set to (0, provided value).
#'
#' #' @return List of pulsar result objects as produced by method analyse_pulsar_data. Each contains
#' plot, peak descriptors, peak summary and animal summary for one data series.
#'
#' @export
#===================================================================================
run_pulsar_batch <- function(input_folder_name,
                             param_list = NULL,
                             make_outfiles = TRUE,
                             y_axis_max = NULL)
{
  # User can set up a param_list if they wish. Otherwise use last stored value,
  # or defaults, in no stored Rdata object
  if (is.null(param_list)) {
    param_list <- restore_param_list() # Last used or default
  }

  # !!! Change to passing list to pulsar method?
  smoothing_fraction <- param_list[["smoothing_fraction"]]
  g_values <- param_list[["g_values"]]
  extinction_threshold <- param_list[["extinction_threshold"]]
  peak_split_depth <- param_list[["peak_split_depth"]]
  sdr_coef = param_list[["sdr_coef"]]
  nearest_nadir_distance = param_list[["nearest_nadir_distance"]]

  # !!! These currently not in use while smoothing algorithm is being discussed
  n_steps = param_list[["n_steps"]]
  if (is.null(n_steps)){
    n_steps = 6
  }
  sdr_assay <- param_list[["sdr_assay"]]

  # Combine all files in folder into a correctly formatted df and experiment name (in list)
  input_data <- read_pulsar_input_folder(input_folder_name)

  # If something goes wrong, read_pulsar_input_folder will return null
  if (is.null(input_data)){
    warning(paste("Error processing data files in", input_folder_name))
    return(NULL)
  }

  # Process as normal
  experiment_name <- input_data$ExperimentName
  hormone_df <- input_data$HormoneDF

  all_results <- analyse_pulsar_data(hormone_df,
                                     smoothing_fraction,
                                     g_values,
                                     extinction_threshold,
                                     peak_split_depth,
                                     sdr_coef,
                                     nearest_nadir_distance,
                                     sdr_assay,
                                     n_steps,
                                     y_axis_max)

  # Write individual peaks and summary peaks to file and generates png for each plot
  if (make_outfiles) {
    write_multi_pulsar(experiment_name, param_list, all_results)
  }


  return (all_results)
} # end run_pulsar_batch





