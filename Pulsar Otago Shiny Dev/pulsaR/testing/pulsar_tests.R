# Clean up and reload
rm(list = ls())

#library(devtools)
#devtools::load_all()

# Source the R script files in the local package copy
package_path <- "pulsaR/R/"
source_files <- c("pulsar_constants.R", "pulsar_utilities.R", "pulsar_computation.R", "pulsar_main.R")
for (source_file in source_files)
{
  source(paste(package_path, source_file, sep=""))
}




#=====================================================================
# Split peak bug fix testing. Deep Debug


data_file_path = "test_data/Five Peak Test.csv"
input_data <- read_pulsar_input_file(data_file_path)

experiment_name <- input_data$ExperimentName
hormone_df <- input_data$HormoneDF

smoothing_fraction <- 0.70
g_values <- c(6.2, 4.9, 2.5, 1.5, 1.2)
extinction_threshold <- 0.1
peak_split_depth <- 1.3
sdr_coef <- c(0, 2.5, 3.3)
nearest_nadir_distance <- 3
n_steps = const$n_steps_constant
                              

check_data_integrity(hormone_df)
n_animals <- ncol(hormone_df) - const$n_utility_columns
sample_times <- as.numeric(hormone_df[,2])
  
# Prepare loop driver sequence from 3 to n, the animal columns
animal_column_indices <- (const$n_utility_columns + 1):(const$n_utility_columns + n_animals)
a <- animal_column_indices[1]
  
animal_id <- colnames(hormone_df)[a]
animal_conc <- as.numeric(hormone_df[,a])
animal_df <- data.frame(time = sample_times, conc = animal_conc)

#####################################################
# Entering code for function pulsar_one
#####################################################

colnames(animal_df) = c("Time", "Concentration")
animal_df <- clean_data(animal_df, extinction_threshold)

times <- animal_df$Time
raw_concentrations <- animal_df$Concentration
smoothed_concentrations <- gen_smoothed(times, raw_concentrations, smoothing_fraction, n_steps)
rescaled_conc <- gen_rescaled(raw_concentrations, smoothed_concentrations, sdr_coef = sdr_coef)
peak_points <- find_peaks(rescaled_conc, g_values)

################################################################################
# Here's the trouble spot...

split_peak_points <- fsm_split_peaks(peak_points, peak_split_depth, rescaled_conc)

################################################################################


trimmed_peak_points <- trim_end_peaks(split_peak_points, rescaled_conc)

debug_df <- animal_df
debug_df$smoothed <- smoothed_concentrations
debug_df$rescaled <- rescaled_conc
debug_df$raw_pulse_flag <- peak_points$flag
debug_df$raw_pulse_length <- peak_points$score
debug_df$split_pulse_flag <- split_peak_points$flag
debug_df$split_pulse_length <- split_peak_points$score
debug_df$trimmed_pulse_flag <- trimmed_peak_points$flag


filename <- paste0(substr(data_file_path,11,25), ".csv")
write.csv(debug_df, filename)



peak_features <- gen_peak_features(trimmed_peak_points, rescaled_conc)

# If you found any peaks, process them. Otherwise leave these elements NULL
peak_descriptors <- NULL
peak_summary <- NULL

gen_peak_descriptors(animal_id,
                                             times,
                                             raw_concentrations,
                                             smoothed_concentrations,
                                             peak_features,
                                             nearest_nadir_distance)

gen_peak_plot(animal_id,
                             animal_df,
                             smoothed_concentrations,
                             trimmed_peak_points,
                             peak_features,
                             nearest_nadir_distance)



































# Test file read in
input_file <- "testing/test_data/orig_data.csv"
read_pulsar_input_file(input_file)


# Test run pulsar on default params (remove param_list.Rdata if present)
# All points will be in peaks
unlink("param_list.Rdata")
input_file_name <- "testing/test_data/acth_data.csv"
param_list = NULL
make_outfiles = FALSE
zip_files = FALSE
y_axis_max = NULL

run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)

# Test on provided parameters
# Only first two high points counted as peaks with high sdr_coef
save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(4, 6, 10),
                nearest_nadir_distance = 3)

run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)

# Test restore param_list. Should be values above
p <- restore_param_list()
p

# Test restore without saved data object. Should be 0/1's from constants
unlink("param_list.Rdata")
p <- restore_param_list()
p



# Test numeric match against deployed shinyapps.io version
# Test input file with multiple animal columns
input_file_name <- "testing/test_data/pulsar_sim_data_06_animals.csv"
param_list = NULL
make_outfiles = FALSE
zip_files = FALSE
y_axis_max = NULL

save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0.4,
                peak_split_depth = 2.7,
                sdr_coef = c(4, 5, 6),
                nearest_nadir_distance = 3)

all_results <- run_pulsar(input_file_name = input_file_name,
                          param_list = param_list,
                          make_outfiles = make_outfiles,
                          zip_files = zip_files,
                          y_axis_max = y_axis_max)


# Animal_ID 	Peak_Time 	Max_Conc 	Amplitude
# arszm 	5.00 	1.81 	0.55
# arszm 	80.00 	4.86 	4.03
# arszm 	185.00 	4.37 	3.91

arszm_id <- 6
all_results[[arszm_id]]$peak_descriptors


# Test writing output files
input_file_name <- "testing/test_data/acth_data.csv"
param_list = NULL
make_outfiles = TRUE
zip_files = FALSE
y_axis_max = NULL

save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(4, 6, 10),
                nearest_nadir_distance = 3)

run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)



# Test input file with multiple animal columns
input_file_name <- "testing/test_data/pulsar_sim_data_06_animals.csv"
param_list = NULL
make_outfiles = FALSE
zip_files = FALSE
y_axis_max = NULL

save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(4, 6, 10),
                nearest_nadir_distance = 3)

all_results <- run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)


# No peaks
check_animal <- 2
print(all_results[[check_animal]]$peak_plot)
print(all_results[[check_animal]]$peak_descriptors)
print(all_results[[check_animal]]$peak_summary)
print(all_results[[check_animal]]$animal_summary)

# Two peaks
check_animal <- 5
print(all_results[[check_animal]]$peak_plot)
print(all_results[[check_animal]]$peak_descriptors)
print(all_results[[check_animal]]$peak_summary)
print(all_results[[check_animal]]$animal_summary)





# Test 3 - batch (change make_outfiles to test)
save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(4,5,6),
                nearest_nadir_distance = 3)

input_folder_name <- "testing/test_data/Batch_Demo"
param_list = NULL
make_outfiles = FALSE
y_axis_max = NULL

all_results <- run_pulsar_batch(input_folder_name = input_folder_name,
                                param_list = param_list,
                                make_outfiles = make_outfiles,
                                y_axis_max = y_axis_max)


check_animal <- 3
print(all_results[[check_animal]]$peak_plot)
print(all_results[[check_animal]]$peak_descriptors)
print(all_results[[check_animal]]$peak_summary)
print(all_results[[check_animal]]$animal_summary)

check_animal <- 5
print(all_results[[check_animal]]$peak_plot)
print(all_results[[check_animal]]$peak_descriptors)
print(all_results[[check_animal]]$peak_summary)
print(all_results[[check_animal]]$animal_summary)


# Adjusted axis
save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(1,1,2),
                nearest_nadir_distance = 3)

input_folder_name <- "testing/test_data/Axis_Tests"
y_axis_max = NULL
all_results <- run_pulsar_batch(input_folder_name = input_folder_name,
                                param_list = param_list,
                                make_outfiles = make_outfiles,
                                y_axis_max = y_axis_max)
all_plots <- c()
for (i in 1:4)
{
  all_plots[[i]] <- all_results[[i]]$peak_plot
}
gridExtra::grid.arrange(grobs = all_plots, nrow = 2)



# Fixed axis
save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(1,1,2),
                nearest_nadir_distance = 3)

input_folder_name <- "testing/test_data/Axis_Tests"
y_axis_max = 50

all_results <- run_pulsar_batch(input_folder_name = input_folder_name,
                                param_list = param_list,
                                make_outfiles = make_outfiles,
                                y_axis_max = y_axis_max)
all_plots <- c()
for (i in 1:4)
{
  all_plots[[i]] <- all_results[[i]]$peak_plot
}
gridExtra::grid.arrange(grobs = all_plots, nrow = 2)


# Same data and params as Merriam
input_file_name <- "testing/test_data/orig_data.csv"
param_list = NULL
make_outfiles = FALSE
zip_files = FALSE
y_axis_max = NULL

save_parameters(smoothing_fraction = 1,
                g_values = c(5, 2.6, 1.9, 1.5, 1.2),
                extinction_threshold = 0.4,
                peak_split_depth = 2.7,
                sdr_coef = c(0,0,8),
                nearest_nadir_distance = 6)

orig_results <- run_pulsar(input_file_name = input_file_name,
                           param_list = param_list,
                           make_outfiles = make_outfiles,
                           zip_files = zip_files,
                           y_axis_max = y_axis_max)

print(orig_results[[1]]$peak_plot)


# Trimming ambiguous points
save_parameters(smoothing_fraction = 0.75,
                g_values = c(3, 2.5, 2.0, 1.5, 1.0),
                extinction_threshold = 0,
                peak_split_depth = 2.7,
                sdr_coef = c(4,7,9),
                nearest_nadir_distance = 3)

input_folder_name <- "testing/test_data/Trim_Tests"
param_list = NULL
make_outfiles = FALSE
y_axis_max = NULL

all_results <- run_pulsar_batch(input_folder_name = input_folder_name,
                                param_list = param_list,
                                make_outfiles = make_outfiles,
                                y_axis_max = y_axis_max)

# Should have one ambiguous and two real peaks
check_animal <- 1
print(all_results[[check_animal]]$peak_plot)

print(all_results[[check_animal]]$peak_descriptors)
print(all_results[[check_animal]]$peak_summary)
print(all_results[[check_animal]]$animal_summary)

all_plots <- c()
for (i in 1:4)
{
  all_plots[[i]] <- all_results[[i]]$peak_plot
}
gridExtra::grid.arrange(grobs = all_plots, nrow = 2)


# Test writing to zip file for remote download
# input_file <- "testing/test_data/pulsar_sim_data_06_animals.csv"
# save_parameters(smoothing_fraction = 0.75,
#                 g_values = c(2.5, 2.0, 1.9, 1.7, 1.6),
#                 extinction_threshold = 0.4,
#                 peak_split_depth = 2.7,
#                 sdr_coef = c(4,7,9),
#                 nearest_nadir_distance = 3)
# run_pulsar(input_file, make_outfiles = TRUE, zip_files = TRUE)





# Test near-nadir amplitude computation
input_file_name <- "testing/test_data/nadir_amplitude_01.csv"
param_list = NULL
make_outfiles = FALSE
zip_files = FALSE
y_axis_max = NULL


save_parameters(smoothing_fraction = 0.75,
                g_values = c(2.5, 2.0, 1.9, 1.7, 1.6),
                extinction_threshold = 0.4,
                peak_split_depth = 2.7,
                sdr_coef = c(4,7,9),
                nearest_nadir_distance = 1)

run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)

save_parameters(smoothing_fraction = 0.75,
                g_values = c(2.5, 2.0, 1.9, 1.7, 1.6),
                extinction_threshold = 0.4,
                peak_split_depth = 2.7,
                sdr_coef = c(4,7,9),
                nearest_nadir_distance = 3)

run_pulsar(input_file_name = input_file_name,
           param_list = param_list,
           make_outfiles = make_outfiles,
           zip_files = zip_files,
           y_axis_max = y_axis_max)



