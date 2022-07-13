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
#' read_pulsar_input_file
#'
#' Read and prepare a single pulsaR input data file (see format requirements below)
#'
#' Input file must be .csv.
#'
#' File format is:
#'
#' \itemize{
#' \item Line 1: Experiment name (any string) as the first element
#' \item Line 2: Column header line. Column headers can be any strings. Data
#' series column headers are treated as animal IDs, so should be unique.
#' \item Lines 3 to 500: Data values
#'}
#' Columns should be:
#' \itemize{
#' \item Column 1: Sample numbers (1 to n)
#' \item Column 2: Times (in minutes)
#' \item Columns 3 to n: Data values (typically hormone readings). Each data series (typically
#' animal) is in a single column.
#'}
#' All data series in a single file use the same time values.

#'
#' @param  pulsar_file_name Complete file path to a correctly formatted input
#'   file
#'
#' @return A named list (ExperimentName, HormoneDF) containing the experiment
#'   name and a data frame in correct format for further pulsaR processing. Returns
#'   NULL on I/O or data formatting error.
#'
#' @import utils
#' @export
#================================================================================

read_pulsar_input_file <- function(pulsar_file_name)
{
  # Prepare for check condition on read fail
  hormone_df <- NULL

  # Drop out on fail to keep the conditional structure simpler
  if (!file.exists(pulsar_file_name)) {
   return(NULL)
  }

  # R's awkward trycatch syntax...
  format_check <- tryCatch(
    # Read in the csv file, suppressing R's desire to put an X in front of colnames it doesn't like
    {  hormone_df <- utils::read.csv(pulsar_file_name, check.names = FALSE, )},
    error = function(e) return (NULL)
  )

  if (is.null(format_check)){
    return(NULL)
  }


  # If the file is empty, or otherwise returns Null on read,
  # drop out
  if (is.null(hormone_df)) {
    return (NULL)
  }

  # Check the structure of the returned csv here, to avoid crashes

  drop_out <- FALSE
  # Must have at least 3 rows and at least 3 columns to be viable
  if ((ncol(hormone_df) < 3) | (nrow(hormone_df) < 3)){
    drop_out <- TRUE
  }

  # !!! TODO: Add other checks here....
  if (drop_out){
    return(NULL)
  }

  # If we reach here, it looks ok to try to parse the file.
  # Data integrity checks can be performed outside this method.

  # The header line in the file becomes the column names of the df
    # That's how read.csv works
    header_line <- colnames(hormone_df)[1]

    # The column names from the file become the first row of the df
    # That's a requirement of the Pulsar input file format
    orig_col_names <- hormone_df[1,]

    # Delete that first row
    hormone_df <- hormone_df[-1,]

    # Cast all readings to numeric (read.csv imports as chr if col contains chr)
    hormone_df <- data.frame(lapply(hormone_df, as.numeric))

    # Reassign the column names
    colnames(hormone_df) <- orig_col_names

    # Return the prepared df with the experiment name and the data df
    return(list(ExperimentName = header_line, HormoneDF = hormone_df))

} # end read_pulsar_input_file


#================================================================================
#' read_pulsar_input_folder
#'
#' Read and prepare multiple pulsaR input data files (each containing exactly
#' one data series; see format requirements below) contained in a single folder.
#'
#' Input files must be .csv
#'
#' File format is:
#'
#' \itemize{
#' \item Line 1: Experiment name (any string) as the first element
#' \item Line 2: Column header line. Data series column header is treated as animal ID
#' so should be unique across all files in the folder.
#' \item Lines 3 to 500: Data values
#'}
#' Columns should be:
#' \itemize{
#' \item Column 1: Sample numbers (1 to n)
#' \item Column 2: Times (in minutes)
#' \item Column 3 Data values (typically hormone readings) for a single data series (typically
#' an animal)
#'}
#' All data files in the folder must have the same time values.

#'
#' @param  input_folder_name Complete file path to a correctly formatted input
#'   file
#'
#' @return A named list (ExperimentName, HormoneDF) containing the experiment
#'   name and a data frame combining all files in correct format for further
#'   pulsaR processing. Returns NULL on I/O or data formatting error
#'
#' @import utils
#' @export
#================================================================================
read_pulsar_input_folder <- function(input_folder_name)
{


    all_input_files = NULL

    # If the input arg is null...
    if (is.null(input_folder_name)){
      warning("Please supply an input folder name")
      return(NULL)
    }

    # If the folder can't be found...
    if (!file.exists(input_folder_name)){
      warning("Input folder not found")
      return(NULL)
    }

    # Try to read the folder contents
    all_input_files <- list.files(input_folder_name)


    # User must provide either an available file name or a folder name
    if (is.null(all_input_files) | length(all_input_files) == 0) {
      warning("Folder contains no input files")
      return(NULL)
    }

    # All times must be the same, so we can get them from any of the files
    # We use the first one
    full_first_file_name <- paste(input_folder_name, all_input_files[1], sep="/")
    first_file <- read_pulsar_input_file(full_first_file_name)
    first_df <- first_file$HormoneDF

    # Same for the sample numbers -- must be identical for all files in the folder
    sample_numbers <- as.numeric(first_df[[1]])
    times <- as.numeric(first_df[[2]])

    # Prepare data frame to accrue
    hormone_df <- data.frame(Sample = sample_numbers,
                             Time = times)


    # Gather all files into df

    # !!! Currently batch files contain only one data column. Will start in column 3 after Sample and Time
    data_column_index <- 3

    # Keep track so we can assign the column headers.
    animal_column_count <- data_column_index

    for (input_file_name in all_input_files)
    {

      full_file_name <- paste(input_folder_name, input_file_name, sep="/")

      # Returns list with ExperimentName and HormoneDF
      input_data <- read_pulsar_input_file(full_file_name)

      # If there was a problem trying to read in the data file,
      # read_pulsar_input_file returns NULL
      if (is.null(input_data)){
        warning(paste("File", full_file_name, "did not import correctly"))
        return(NULL)
      }

      # Grab the concentrations column
      input_df <- input_data$HormoneDF
      input_concentrations <- input_df[[data_column_index]]

      # All data columns must be the same length as the times column
      if (length(input_concentrations) != length(times)){
        warning(paste("File", full_file_name, ": data column", data_column_index, "has the wrong number of values"))
        return(NULL)
      }


      # Add it to the growing df
      hormone_df$newColumn <- as.numeric(input_concentrations)

      # Change the column name to the animal id. In batch mode, this is
      # the input file name with its file suffix stripped off
      animal_id_col_name <- gsub(".csv", "", input_file_name)
      colnames(hormone_df)[animal_column_count] <- animal_id_col_name

      # Ready for the next one
      animal_column_count <- animal_column_count + 1
    } # end for all files in folder


    # Prepare the returned object
    # For batch processing, the experiment name is set to the folder name

    # Strip any path prefix
    experiment_name <- basename(input_folder_name)
    pulsar_input_list <- list(ExperimentName = experiment_name,
                              HormoneDF = hormone_df)
    return(pulsar_input_list)
} # end read folder


#================================================================================
#' clean_data
#'
#' Accepts a df containing a single data series. Removes negative readings.
#' Replaces non-negative readings below extinction threshold with threshold value.
#' Returns the cleaned version.
#'
#' @param raw_data_df The df is expected to have a column Concentration as
#' produced by read_pulsar_input_file and read_pulsar_input_folder.
#'
#' @param threshold A double defining the minimum value in the cleaned
#' data series. Represents the sensitivity limit of the assay. May be 0.
#'
#' @return Copy of original df with cleaned Concentration column
#'
#' @export
#================================================================================
clean_data <- function(raw_data_df, threshold)
{
  # Select only rows that are not NA
  clean_data_df <- raw_data_df[!is.na(raw_data_df$Concentration),]

  # Select only  rows for which the value of Concentration is > 0
  clean_data_df <- clean_data_df[(clean_data_df$Concentration >= 0),]

  # Replace concentrations below threshold
  clean_data_df$Concentration[clean_data_df$Concentration < threshold] <- threshold

  return(clean_data_df)
}


#================================================================================
#' check_data_integrity
#'
#' Accepts a df with Sample No., Time and Data columns as produced by read_pulsar_input_file
#' and read_pulsar_input_foder. Checks for a variety of conditions that might indicate a problem
#' with the data and issues warnings.
#'
#' @param hormone_df A df with Sample No., Time and Data columns as produced by read_pulsar_input_file
#' and read_pulsar_input_foder.
#'
#' @export
#================================================================================
check_data_integrity <- function(hormone_df)
{

  # !!! TODO: Can add additional checks useful here as they
  # are identified.

  # Hormone_df should have sample number in column 1, times in column 2 and
  # data in columns 3 to n

  sample_numbers <- hormone_df[,1]
  times <- hormone_df[,2]
  readings <- hormone_df[,3:ncol(hormone_df)]

  # Check for suspicious things in sample_numbers. Should be 1 to n.
  sample_template <- 1:length(sample_numbers)
  sample_check <- sample_numbers == sample_template
  if (!(all(sample_check))){
    warning("Sample numbers may be in wrong format")
  }


  # Check for suspicious things in time. Should be monotonically
  # increasing and cannot be negative
  for (i in 1:(length(times)-1))
  {
    if (times[i] >= times[i+1] ){
      warning("Sample times are not monotonically increasing")
    }
  }

  if (any(times < 0)){
    warning("Sample times contains a negative value")
  }

  # Check for suspicious things in the data columns
  # Non-numeric values !!! Doesn't work for one-column
  #if (any(!apply(readings, 2, is.numeric))) {
   # warning("Concentration readings contain non-numeric values")
  #}

  # Negative values; skipping NA
  if (any(readings < 0, na.rm = TRUE)){
    warning("Concentration readings contain negative values")
  }
  # Missing values -- jagged data reads in with NA to fill
  if (any(is.na(readings))){
    warning("Concentration data series has missing values.")

  }


}


#================================================================================
#' get_current_time_string
#'
#' Generates a string time stamp. Used for outfile naming.
#'
#' @export
#================================================================================
get_current_time_string <- function()
{
  now <- as.POSIXlt(Sys.time())

  year <- now$year - 100
  month <- now$mon
  day <- now$mday
  hour <- now$hour
  min <- now$min
  sec <- trunc(now$sec)

  string_time <- paste(year, month, day, hour, min, sec, sep="_")

  return(string_time)
}



#================================================================================
#' save_parameters
#'
#' Saves out a set of the 5 adjustable algorithm parameters as an Rdata
#' object for subsequent reload.
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
#' points before the max point of the peak/pulse event. The nearest_nadir_distance
#' parameter sets the value of n. This value will default to 3, a sensible distance for
#' LH data.
#'
#' @export
#================================================================================
save_parameters <- function(smoothing_fraction,
                            g_values,
                            extinction_threshold,
                            peak_split_depth,
                            sdr_coef,
                            nearest_nadir_distance)
{

  # !!! These currently not adjustable.
  sdr_assay = NULL
  n_steps = const$n_steps_constant

  param_list = list(smoothing_fraction = smoothing_fraction,
                    g_values = g_values,
                    extinction_threshold = extinction_threshold,
                    peak_split_depth = peak_split_depth,
                    sdr_coef = sdr_coef,
                    nearest_nadir_distance = nearest_nadir_distance,
                    sdr_assay = sdr_assay,
                    n_steps = n_steps)

  config_file_name = "param_list.Rdata"
  save(param_list, file = config_file_name)
}

#================================================================================
#' restore_param_list
#'
#' Loads file param_list.Rdata if present, creating object
#'  param_list. If the Rdata file is not present, default values (following Merriam 1982)
#'  are supplied. Later routines can parse param_list for the individual parameter terms.
#'
#'
#' @return Named list containing the six adjustable algorithm parameters (smoothing fraction,
#' g(n) values, extinction threshold, peak split depth, assay variability estimation coefficients
#' and nearest nadir distance), and the two fixed algorithm parameters (constant assay variability
#' term and number of steps for lowess).
#'
#' @export
#================================================================================
restore_param_list <- function()
{
  if (file.exists("param_list.Rdata")) {
    load("param_list.Rdata")
  } else {
    # Default values based on Otago testing. See pulsar_constants.R
    smoothing_fraction <- const$smoothing_fraction
    g_values <- const$g_values
    extinction_threshold <- const$extinction_threshold
    peak_split_depth <- const$peak_split_depth
    sdr_coef = const$sdr_coef
    nearest_nadir_distance = const$nearest_nadir_distance

    sdr_assay = const$sdr_assay # Not used in this implementation
    n_steps = const$n_steps # improves lowess accuracy over original value of 3


    param_list = list(smoothing_fraction = smoothing_fraction,
                      g_values = g_values,
                      extinction_threshold = extinction_threshold,
                      peak_split_depth = peak_split_depth,
                      sdr_coef = sdr_coef,
                      nearest_nadir_distance = nearest_nadir_distance,
                      sdr_assay = sdr_assay,
                      n_steps = n_steps)
  }
  return(param_list)
}


#================================================================================
#' gen_base_pulsar_plot
#'
#' Generates a ggplot of raw concentrations and smoothed concentrations.
#' If there are identified peaks, those points can be added later as
#' additional ggplot layers (see summarise_pulsar_results).
#'
#' @param animal_df Data frame containing three columns: Sample number (ordinal),
#' sample time (in minutes), and raw data values.
#'
#' @param animal_id Unique string identifier for current data series. Used in
#' image title
#'
#' @param smoothed_concentrations Vector containing smoothed (baseline) data values. Must
#' be the same length as the raw data column in animal_df. See method gen_smoothed
#'
#' @param time_column_label Axis label for x axis
#'
#' @param concentration_column_label Axis label for y axis
#'
#' @param pulsar_font Font object to insure good cross-platform behaviour. See
#' file pulsar_constants.R.
#'
#' @param y_axis_max Allows the user to fix the upper limit of the y-axis. If NULL
#' y-axis is adjusted to the range of the raw concentration data values.
#'
#'
#' @export
#================================================================================
gen_base_pulsar_plot <- function(animal_df, animal_id, smoothed_concentrations, time_column_label,
                                 concentration_column_label, pulsar_font, y_axis_max = NULL)
{
    peak_plot <- ggplot2::ggplot(data = animal_df, ggplot2::aes_string(x = time_column_label)) +
      ggplot2::geom_point(ggplot2::aes_string(y = concentration_column_label)) +
      ggplot2::geom_line(ggplot2::aes_string(y = concentration_column_label), color = "black") +
      ggplot2::geom_line(ggplot2::aes(y = smoothed_concentrations), color = "blue") +
      ggplot2::ggtitle(animal_id) +
      ggplot2::theme_bw() +
      ggplot2::theme(text=ggplot2::element_text(family=pulsar_font, size = 14))

    if (!is.null(y_axis_max))
    {
      peak_plot <- peak_plot +
        ggplot2::ylim(0, y_axis_max)
    }

    return(peak_plot)
}


#================================================================================
#' write_multi_pulsar
#'
#' Generates output files from a complete pulsar analysis.
#' Accepts a list of pulsar results ($plot, $peak_descriptors & $peak_summary) as produced
#' by end-user-facing function run_pulsar and run_pulsar_batch.
#' Generates combined output files for all data series and individual images files
#' for each data series.
#'
#'
#' @param  experiment_name String value for experiment name. Typically one
#' uses the initialised value from reading and input file or folder.
#' @param param_list Named list containing algorithm parameters as produced by
#' restore_param_list
#' @param pulsar_result_list A list of pulsar results  (named fields $plot,
#' $peak_descriptors & $peak_summary) as produced by functions run_pulsar and
#' run_pulsar_batch.
#' @param zip_files Set to true to generate single zip containing all outfile png and csv.
#' Can be used for remote download when deployed as a shiny app.
#'
#' @import utils
#' @import gridExtra
#' @import grDevices
#' @export
#================================================================================
write_multi_pulsar <- function(experiment_name, param_list, pulsar_result_list, zip_files = FALSE)
{

  # Make folders to hold output files
  # If downloading as a zip (i.e. from a remote web host) just write everything to the working directory
  # Otherwise set up separate folders for the images and the csv files
  if (zip_files) {
    image_folder_name <- getwd();
    analysis_folder_name <- getwd();
  } else {
    # For images
    image_folder_name <- paste(experiment_name, "Pulsar Images")
    if (!file.exists(image_folder_name)) {
     dir.create(image_folder_name)
    }
    # For results files
    analysis_folder_name <- paste(experiment_name, "Pulsar Analysis")
    if (!file.exists(analysis_folder_name)) {
     dir.create(analysis_folder_name)
    }
  }


  # Set up dfs for hold tables for output
  # indiv_peaks_df <- data.frame(AnimalID = list(),
  #                              AvgStartTime = list(),
  #                              PeakTime = list(),
  #                              MaxConc = list(),
  #                              Amplitude = list(),
  #                              PeakLength = list())

  # 19-06-21 Start time and Peak Length removed as client doesn't want them
  indiv_peaks_df <- data.frame(Animal_ID = list(),
                               Peak_Time = list(),
                               Max_Conc = list(),
                               Amplitude = list())


  summary_peaks_df <- data.frame(Animal_ID = list(),
                                 Avg_Amp = list(),
                                 SD_Amp = list(),
                                 Avg_Peak_Interval = list(),
                                 SD_Peak_Interval = list())

  summary_animals_df <- data.frame(Animal_ID = list(),
                                   Avg_Concentration = list(),
                                   SD_Concentration = list())



  n_results <- length(pulsar_result_list)

  current_time_string <- get_current_time_string()


  # Save the parameter values to a readable file
  # Awkward because the list is jagged

  serialised_parameters <- data.frame(smoothing_fraction = param_list$smoothing_fraction,
                                      g1 = param_list$g_values[1],
                                      g2 = param_list$g_values[2],
                                      g3 = param_list$g_values[3],
                                      g4 = param_list$g_values[4],
                                      g5 = param_list$g_values[5],
                                      extinction_threshold = param_list$extinction_threshold,
                                      peak_split_depth = param_list$peak_split_depth,
                                      sdrQuad = param_list$sdr_coef[1],
                                      sdrLin = param_list$sdr_coef[2],
                                      sdrConst = param_list$sdr_coef[3],
                                      nearest_nadir_distance = param_list$nearest_nadir_distance)

  parameter_file_name <- paste(analysis_folder_name, "/",
                               experiment_name, "_parameters_", current_time_string, ".csv", sep="")

  utils::write.csv(serialised_parameters, parameter_file_name)


  # Walk through all the results accruing contents into the data frames and generating the
  # individual images
  for (r in 1:n_results)
  {
    animal_output <- pulsar_result_list[[r]]
    animal_id <- animal_output$animal_summary$Animal_ID

    indiv_peaks_df <- rbind(indiv_peaks_df, animal_output$peak_descriptors)
    summary_peaks_df <- rbind(summary_peaks_df, animal_output$peak_summary)
    summary_animals_df <- rbind(summary_animals_df, animal_output$animal_summary)

    image_file_name <- paste(image_folder_name, "/",
                             experiment_name, "_", animal_id, "_", current_time_string, ".png", sep="")

    grDevices::png(image_file_name)
    print(animal_output$peak_plot)
    grDevices::dev.off()
  }

  # Making a single image with all the plots when there are multiple

  if (n_results > 1) {
    all_plots <- list()
    for (r in 1:n_results)
    {
      animal_output <- pulsar_result_list[[r]]
      all_plots[[r]] <- animal_output$peak_plot
    }

    all_image_file_name <- paste(image_folder_name, "/",
                                 experiment_name, "_all_", current_time_string, ".png", sep="")

    grDevices::png(all_image_file_name, width = const$mult_image_width, height = const$mult_image_height)

    ## Use grid.arrange for ggplots. par doesn't work
    gridExtra::grid.arrange(grobs = all_plots, nrow = 2)
    grDevices::dev.off()
  }

  # Save out the spreadsheets
  peak_ind_outfile <- paste(analysis_folder_name, "/",
                            experiment_name, "_peak_desc_", current_time_string, ".csv", sep="")
  peak_sum_outfile <- paste(analysis_folder_name, "/",
                            experiment_name, "_peak_summ_", current_time_string, ".csv", sep="")
  animal_sum_outfile <- paste(analysis_folder_name, "/",
                              experiment_name, "_conc_summ_", current_time_string, ".csv", sep="")

  utils::write.csv(indiv_peaks_df, peak_ind_outfile)
  utils::write.csv(summary_peaks_df, peak_sum_outfile)
  utils::write.csv(summary_animals_df, animal_sum_outfile)

  # Gather up all the files we wrote. We'll use this to zip, if requested,
  # and in either case return from the method in case caller needs to know.
  all_files <- list.files(path = getwd(), pattern = ".csv$")
  all_files <- c(all_files, list.files(path = getwd(), pattern = ".png$"))

  # If zipping, handle here
  if (zip_files)
  {
    zip_file_name <- paste(experiment_name, ".zip", sep="")
    zip_status <- zip::zip(zip_file_name, all_files)

    # !!! Deal with status == 127, which is the error code
  }

  return(all_files)

}
