# Clean up and reload
rm(list = ls())
library(devtools)
devtools::load_all()


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



