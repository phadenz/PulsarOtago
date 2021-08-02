# https://phaden.shinyapps.io/pulsaR_Otago/

#================================================================================
#     R implementation of PULSAR algorithm (Merriam & Wachter, 1982)

#     Copyright (C) 2021, University of Otago, Patricia Haden, RTIS
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



#=========================================================================
# Set up
#=========================================================================

# Helper function for ui
inline_numericInput=function(ni){
  tags$div( class="form-inline",ni)
}

options(warn = -1)

# Install packages if necessary
if (!require(ps)) install.packages('ps')
library(ps)

if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

if (!require(shinyFeedback)) install.packages('shinyFeedback')
library(shinyFeedback)

if (!require(shinyjs)) install.packages('shinyjs')
library(shinyjs)

if (!require(shinyBS)) install.packages('shinyBS')
library(shinyBS)

if (!require(shinyFiles)) install.packages('shinyFiles')
library(shinyFiles)

# tooltip_text.R contains the string constants
source("tooltip_text.R")

# Source the R script files in the local package copy
package_path <- "pulsaR/R/"
source_files <- c("pulsar_constants.R", "pulsar_utilities.R", "pulsar_computation.R", "pulsar_main.R")
for (source_file in source_files)
{
  source(paste(package_path, source_file, sep=""))
}

# Layout constants
button_width <- "250px"
g_box_width <- "150px"



#=========================================================================
# Utility functions
#=========================================================================

#=======================================
# Update input controls
update_input_controls <- function(session)
{
  # On launch, set up params and values if running locally. 
  # (When running remotely, just keep defaults). 
  # Method restore_param_list loads param_list.Rdata, if it's
  # there, otherwise assigns defaults.
  param_list <- restore_param_list()
  
  # !!! Change to passing list to pulsar method?
  smoothing_fraction <- param_list[["smoothing_fraction"]]
  g_values <- param_list[["g_values"]]
  extinction_threshold <- param_list[["extinction_threshold"]]
  peak_split_depth <- param_list[["peak_split_depth"]]
  sdr_coef = param_list[["sdr_coef"]]
  n_steps = param_list[["n_steps"]]
  sdr_assay <- param_list[["sdr_assay"]]
  
  updateTextInput(session, "spn_smoothing_fraction", value = smoothing_fraction)
  updateTextInput(session, "spn_split_depth", value = peak_split_depth)
  updateTextInput(session, "spn_level_detection", value = extinction_threshold)
  updateTextInput(session, "spn_g1", value = g_values[1])
  updateTextInput(session, "spn_g2", value = g_values[2])
  updateTextInput(session, "spn_g3", value = g_values[3])
  updateTextInput(session, "spn_g4", value = g_values[4])
  updateTextInput(session, "spn_g5", value = g_values[5])
  updateTextInput(session, "spn_quad", value = sdr_coef[1])
  updateTextInput(session, "spn_lin", value = sdr_coef[2])
  updateTextInput(session, "spn_const", value = sdr_coef[3])
  
  return (param_list)
  
}

#=======================================
clear_tables <- function(output)
{
  # Clear tables
  output$tbl_individual_peaks <- renderTable(NULL)
  output$tbl_peak_summary <- renderTable(NULL)
  output$tbl_animal_summary <- renderTable(NULL)
}

#=======================================
clear_display_controls <- function(input, output, session)
{
  # When they load a new data file, clear everything away
  # !!! Tried shinyjs::show/hide, but the page reflow looks messy
  
  # Clear tabset
  output$placeholder_tbs_images <- renderUI({h4("")})
  
  # Clear plots and tab labels
  null_plot <- NULL
  null_tab_label <- ""
  
  max_animals <- 9
  max_indices <- 1:max_animals
  lapply(max_indices,
         FUN = function(i)
         {
           plot_id <- paste("plt_0", i, sep="")
           output[[plot_id]] <- renderPlot(null_plot)
           #shinyjs::hide(plot_id)
           
           tab_id <- paste("txt_panel_0", i, sep="")
           output[[tab_id]] <- renderText(null_tab_label)
           #shinyjs::hide(tab_id)
         }
  )
  
  clear_tables(output)
  
  
} # end clear_display_controls


#=======================================
render_tooltip_triggers <- function(output, session)
{
  trigger <- "<font color=\"#CA2F0E\"><b>?</b>"
  output$q_smoothing_fraction <- renderText(trigger)
  output$q_split_depth <- renderText(trigger)
  output$q_level_detection <- renderText(trigger)
  output$q_g_function <- renderText(trigger)
  output$q_assay_variability <- renderText(trigger)
  output$q_nearest_nadir_distance <- renderText(trigger)
}

#=======================================
# Allows 9 tabs for easier name parsing
rebuild_tab_panel <- function(input, output, session)
{
  output$placeholder_tbs_images <- renderUI({
    tabsetPanel(id = "tbs_images", 
                tabPanel(title = textOutput("txt_panel_01"), 
                         plotOutput("plt_01"), value = 1),
                tabPanel(title = textOutput("txt_panel_02"), 
                         plotOutput("plt_02"), value = 2),
                tabPanel(title = textOutput("txt_panel_03"), 
                         plotOutput("plt_03"), value = 3),
                tabPanel(title = textOutput("txt_panel_04"), 
                         plotOutput("plt_04"), value = 4),
                tabPanel(title = textOutput("txt_panel_05"), 
                         plotOutput("plt_05"), value = 5),
                tabPanel(title = textOutput("txt_panel_06"), 
                         plotOutput("plt_06"), value = 6),
                tabPanel(title = textOutput("txt_panel_07"), 
                         plotOutput("plt_07"), value = 7),
                tabPanel(title = textOutput("txt_panel_08"), 
                         plotOutput("plt_08"), value = 8),
                tabPanel(title = textOutput("txt_panel_09"), 
                         plotOutput("plt_09"), value = 9)
    )
    
  })
} # end rebuild_tab_panel

#=======================================
clear_and_draw_plots <- function(input, output, session, all_results)
{
  # Clean out any old plots, from earlier runs of the same input file with different
  # arg values !!! NB: for loop doesn't work here. All passes
  # take value from final pass. Shiny weirdness
  
  null_plot <- NULL
  null_tab_label <- ""
  
  max_animals <- 9
  max_indices <- 1:max_animals
  lapply(max_indices,
         FUN = function(i)
         {
           plot_id <- paste("plt_0", i, sep="")
           output[[plot_id]] <- renderPlot(null_plot)
           
           tab_id <- paste("txt_panel_0", i, sep="")
           output[[tab_id]] <- renderText(null_tab_label)
         }
  )
  
  
  
  # Draw the plots on the appropriate tabs, and set the tab labels
  n_animals <- length(all_results)
  animal_indices <- 1:n_animals
  lapply(animal_indices,
         FUN = function(i)
         {
           plot_id <- paste("plt_0", i, sep="")
           tab_id <- paste("txt_panel_0", i, sep="")
           
           if (!is.null(all_results[[i]])) {
             output[[plot_id]] <- renderPlot(all_results[[i]]$peak_plot)
             
             output[[tab_id]] <- renderText(all_results[[i]]$animal_summary$Animal_ID[1])
           }
         }
  )
} # end clean_and_draw_plots


#=========================================
# Hard-coded fields on result object list
fill_tables <- function(input, output, session, result_index, all_results)
{
  pulsar_output <- all_results[[result_index]]
  output$tbl_individual_peaks <- renderTable(pulsar_output$peak_descriptors)
  output$tbl_peak_summary <- renderTable(pulsar_output$peak_summary)
  output$tbl_animal_summary <- renderTable(pulsar_output$animal_summary)
  
  shinyjs::show("tbl_individual_peaks")
  shinyjs::show("tbl_peak_summary")
  shinyjs::show("tbl_animal_summary")
} # end fill_tables



#=========================================================================
# ui
#=========================================================================

ui <- fluidPage(

  shinyjs::useShinyjs(),
  shinyFeedback::useShinyFeedback(),
  
  # Not clear why these have to be separate, but they don't work when bundled in with the head tags....
  tags$style(type='text/css', "#btn_pulsar { width:100%; margin-top: 25px; border: 2px solid #aaaaaa;}"),
  tags$style(type='text/css', "#btn_batch { width:100%; border: 2px solid #aaaaaa;}"),
  tags$style(type='text/css', "#btn_save { width:100%; border: 2px solid #aaaaaa;}"),
  tags$style(type='text/css', "#btn_download { width:100%; border: 2px solid #aaaaaa;}"),
  

  tags$head(tags$style("#side_panel{
                       padding-left:10px;
                       }
                       .form-group {
                       margin-bottom: 15px !important;
                       }
                       .form-inline .form-control {
                       width: 100%
                       }
                       .shiny-notification {
                       color:#ffffff;
                       background-color:#112446;
                       font-size: 200%;
                       width: 200%;
                       position: fixed;
                       top: 90%;
                       left: 70%;
                       }
                       hr{
                       border-top: 2px solid #aaaaaa;}
                       }
                       "
                      )),

  titlePanel(title="Pulsar Otago"),

  sidebarLayout(

    sidebarPanel(width = 6, id = "side_panel",
                 
      fluidRow(
        column(12,  a(href="PULSAR_Otago_Input_Output_Formats_v1.4", "Download File Management Documentation", download=NA, target="_blank"))
      ),
                 
      fluidRow(
        column(12, h4("Pulsar Parameters"))
      ),
      
      fluidRow(
       
        column(3, 
               htmlOutput("q_smoothing_fraction"),
               bsPopover("q_smoothing_fraction", title = "", content = tip_spn_smoothing_fraction),
               numericInput("spn_smoothing_fraction", "Smoothing Fraction", 
                            value = 1.0, min = 0, max = 1, step = 0.01, width = "150px")),
        column(3, 
               htmlOutput("q_split_depth"),
               bsPopover("q_split_depth", title = "", content = tip_spn_split_depth),
               numericInput("spn_split_depth", "Peak Split Depth", 
                            value = 0.0, min = 0, max = 5, step = 0.1, width = "150px")),
        column(3, 
               htmlOutput("q_level_detection"),
               bsPopover("q_level_detection", title = "", content = tip_spn_level_detection),
               numericInput("spn_level_detection", "Level of Detection", 
                            value = 0.0, min = 0, max = 5, step = 0.01, width = "150px")),
        column(3, 
               htmlOutput("q_nearest_nadir_distance"),
               bsPopover("q_nearest_nadir_distance", title = "", content = tip_spn_nearest_nadir_distance),
               numericInput("spn_nearest_nadir_distance", "Amplitude distance", 
                            value = 3, min = 0, max = 50, step = 1, width = "150px"))
      ),
      
      h4(),

      fluidRow(column(6,
                      htmlOutput("q_g_function"),
                      bsPopover("q_g_function", title = "", content = tip_g_function),
                      h4("G Values"))
      ),
      
      fluidRow(
          column(2, inline_numericInput(numericInput("spn_g1", "g(1)", value = 0.0, min=0, max=20, step=0.1))),
          column(2, inline_numericInput(numericInput("spn_g2", "g(2)", value = 0.0, min=0, max=20, step=0.1))),
          column(2, inline_numericInput(numericInput("spn_g3", "g(3)", value = 0.0, min=0, max=20, step=0.1))),
          column(2, inline_numericInput(numericInput("spn_g4", "g(4)", value = 0.0, min=0, max=20, step=0.1))),
          column(2, inline_numericInput(numericInput("spn_g5", "g(5)", value = 0.0, min=0, max=20, step=0.1)))
      ),
      
      h4(),
      
      fluidRow(column(6,
                      htmlOutput("q_assay_variability"),
                      bsPopover("q_assay_variability", title = "", content = tip_assay_variability),
                      h4("Assay Variability"))
      ),
      
      fluidRow(
        column(4, numericInput("spn_quad", "Quadratic", value = 0.0, min=0, max=100, step=0.1, width = "150px")),
        column(4, numericInput("spn_lin", "Linear", value = 0.0, min=0, max=100, step=0.1, width = "150px")),
        column(4, numericInput("spn_const", "Constant", value = 1.0, min=0, max=100, step=0.1, width = "150px"))
      ),
      
      hr(),
      
      h4("Options"),
      
      fluidRow(
        column(6, checkboxInput("chk_fix_y_axis", strong("Fix y-axis to maximum of all data series"), 
                                value = FALSE))
      ),
      
      hr(),
      
      h4("Process Single Input File"),
      
      fluidRow(
        column(7, fileInput("fin_data", "Upload Data File",
                             multiple=FALSE,
                             accept = c(".csv"),
                             placeholder = "Demo - test data/pulsar_sim_data_06_animals.csv")),
        column(5, actionButton("btn_pulsar", "Run Pulsar Current File", width = button_width))
        
      ),
      
      h4("Process Folder (Only available when running locally. First nine series displayed.)"),
      
      fluidRow(
        column(2, shinyDirButton("fin_batch", label = "Select Folder", title = "Select Folder")),
        column(5, textInput("txt_folder", label = NULL, value = "No Folder Selected")),
        column(5, actionButton("btn_batch", "Run Pulsar Batch", width = button_width)),
      ),
      
      hr(),
      
      fluidRow(
        column(6, downloadButton("btn_download", "Save Results", width = button_width)),
        column(8, actionButton("btn_save", "Save Results", width = button_width))
    ),
    
    # end sidebar layout
    ),

    mainPanel(width = 6,
              fluidRow(
                uiOutput("placeholder_tbs_images")
              ),
      h4(),
      
      fluidRow(tableOutput("tbl_individual_peaks")),
      fluidRow(tableOutput("tbl_peak_summary")),
      fluidRow(tableOutput("tbl_animal_summary"))
    ) # end mainPanel
  ), #end sidebarlayout
  hr(),
  h6("Pulsar Otago. Copyright (c) 2021 University of Otago"),
  h6("This program is licensed under GNU 3 (see http://www.gnu.org/licenses/). This program comes with ABSOLUTELY NO WARRANTY.")
) # end ui <- fluidPage



#=========================================================================
# server <- function(...)
#=========================================================================

server <- function(input, output, session) {
  
  is_local <- Sys.getenv('SHINY_PORT') == ""
  
  # check with rsconnect::showLogs(streaming = TRUE, appName = "pulsaR_Otago")
  cat(file=stderr(), "is_local: ", is_local, "\n" )
  
  #================================================ 
  # On launch, set up tooltips
  render_tooltip_triggers(output, session)
  
  #================================================ 
  # On launch, set up params if running locally
  if (is_local){
    param_list <- update_input_controls(session)
  }
  
  #================================================
  # Decide which save button to show
  #   
  if (is_local){
    show("btn_save")
    hide("btn_download")
  } else {
    hide("btn_save")
    show("btn_download")
  }

  #================================================
  # Data structure for storing dynamic data values
  # Default data file is for testing
  v <- reactiveValues(
    data_file_path = "test data/pulsar_sim_data_06_animals.csv",
    data_folder_path = NULL,
    input_data = NULL,
    experiment_name = NULL,
    hormone_df = NULL,
    all_results = NULL
  )
  

  #================================================
  # Folder selection button trigger setup (shinyBS::)
  
  roots <- getVolumes()()
  shinyDirChoose(input,
                 'fin_batch',
                 roots = roots)
  
  #================================================
  # Folder selection button click event handler
  observeEvent(input$fin_batch,
  {
    # Clear any previous displays
    clear_display_controls(input, output, session);
    
    folder_path <- parseDirPath(roots = roots, input$fin_batch)
    display_folder <- basename(folder_path)
    v$data_folder_path = folder_path
    
    updateTextInput(session, "txt_folder", value = display_folder)
  }) # end observeEvent input$fin_batch
  
  #================================================
  # File selection button click event handler
  observeEvent(input$fin_data,
  {
   # Grab file name
   v$data_file_path <- input$fin_data$datapath
           
   # Clear any previous displays
   clear_display_controls(input, output, session);
                 
  }) # end observeEvent input$fin_data
  
  
  
  #================================================
  # Event handler for image tab set.
  # Get the id of the user-selected tabpanel.
  # If you've got data for it, initialise the main output
  # panel with the appropriate data. Otherwise, hide 
  # everything, effectively disabling extra tabs
  observeEvent(input$tbs_images,
  {
    # Get the id of the user-selected tabpanel
    tab_number <- as.numeric(req(input$tbs_images))

    # Initialise display panel widgets
    if (tab_number <= length(v$all_results))
    {
      pulsar_output <- v$all_results[[tab_number]]
      output$tbl_individual_peaks <- renderTable(pulsar_output$peak_descriptors)
      output$tbl_peak_summary <- renderTable(pulsar_output$peak_summary)
      output$tbl_animal_summary <- renderTable(pulsar_output$animal_summary)
      shinyjs::show("tbl_individual_peaks")
      shinyjs::show("tbl_peak_summary")
      shinyjs::show("tbl_animal_summary")
    } else {
      # Don't have to hide the plots, as the whole tabset is regenerated on file
      # open so the plots are empty
      shinyjs::hide("tbl_individual_peaks")
      shinyjs::hide("tbl_peak_summary")
      shinyjs::hide("tbl_animal_summary")
    }
    
  }) # end  observeEvent input$tbs_images
  

  #================================================
  # Event handler Run Pulsar button
  # Process a single file with 1 or more data series columns.
  # Clear old display values, grab parameters from input widgets,
  # Call pulsar computational methods, storing outputs
  # Update display values
  observeEvent(input$btn_pulsar,
  {
    
    # To allow for different numbers of plots (animals) in different data files within
    # the same session, it works best to just rebuild the tab panel completely. At the 
    # moment we are hard-coding a maximum of 9 tabs (animals) for ease of id management.
    rebuild_tab_panel(input, output, session)
    
    # Clear tables. May be redundant, but the timing can get tangled, so this one is for safety
    clear_tables(output)
  
    
    # Get and store file name. This reactive value is assigned in the file input
    # button even handler
    if (!is.null(input$fin_data$datapath)){
      v$data_file_path <- input$fin_data$datapath
    }
  
        
    # Returns NULL if the import fails. Returns a list ExperimentName and HormoneDF
    # if the import succeeds
    input_data <- read_pulsar_input_file(v$data_file_path)
                  
    if (is.null(input_data)){
      # Show notification and skip remainder of method
      showNotification("Error reading input file", type="error")
      } else {
        
        # Load fields of the reactive data holding objects
        v$input_data <- input_data
        v$experiment_name <- v$input_data$ExperimentName
        v$hormone_df <- v$input_data$HormoneDF
        v$all_results = NULL

      # Grab the current algorithm parameters from the relevant controls
      smoothing_fraction <- input$spn_smoothing_fraction
      g_values <- c(input$spn_g1, input$spn_g2, input$spn_g3, input$spn_g4, input$spn_g5)
      extinction_threshold <- input$spn_level_detection
      peak_split_depth <- input$spn_split_depth
      sdr_coef <- c(input$spn_quad, input$spn_lin, input$spn_const)
      nearest_nadir_distance <- input$spn_nearest_nadir_distance
                 
      # If running locally, save these values to be restored next time.
      if(is_local) {
        save_parameters(smoothing_fraction,
                        g_values,
                        extinction_threshold,
                        peak_split_depth,
                        sdr_coef,
                        nearest_nadir_distance)
      }

      
      # Not using the wrapped function run_pulsar here to get finder control over the input
      # parameters. Could be done though.
      
      # Pop up
      notification_id = showNotification("Finding peaks.")
                  
      # Fixing the y-axes to the maximum value in the data columns (from column 3) if the box is checked
      if (input$chk_fix_y_axis){
        max_raw_conc <- max(v$hormone_df[, c(-1, -2)])
        } else {
          max_raw_conc = NULL
        }
                  
                  
                  
      # Call the main pulsar computation routine, passing required args and storing result
      v$all_results = analyse_pulsar_data(v$hormone_df,
                                          smoothing_fraction,
                                          g_values,
                                          extinction_threshold,
                                          peak_split_depth,
                                          sdr_coef,
                                          nearest_nadir_distance,
                                          y_axis_max = max_raw_conc)
                  
                  
      
      clear_and_draw_plots(input, output, session, v$all_results)            
      
      # First tab is displayed by default. Load up here.
      # Other tabs load on tab click event
      result_index <- 1
      fill_tables(input, output, session, result_index, v$all_results)

      removeNotification(notification_id)

    } # end in input_data not null
  }) # end observeEvent btn_pulsar
  
  
  #=====================================================================
  # Event handler btn_batch
  # Process a folder with multiple single-column input files.
  # !!! TBD Code duplication with single button
  observeEvent(input$btn_batch,
  {
         
   # To allow for different numbers of plots (animals) in different data files within
   # the same session, it works best to just rebuild the tab panel completely. At the 
   # moment we are hard-coding a maximum of 9 tabs (animals) for ease of id management.
   rebuild_tab_panel(input, output, session)
                 
   # Clear tables. May be redundant, but the timing can get tangled, so this one is for safety
   clear_tables(output)
                 
                 
   # Grab the current algorithm parameters from the relevant controls
   smoothing_fraction <- input$spn_smoothing_fraction
   g_values <- c(input$spn_g1, input$spn_g2, input$spn_g3, input$spn_g4, input$spn_g5)
   extinction_threshold <- input$spn_level_detection
   peak_split_depth <- input$spn_split_depth
   sdr_coef <- c(input$spn_quad, input$spn_lin, input$spn_const)
   nearest_nadir_distance <- input$spn_nearest_nadir_distance
                 
   # If running locally, save these values to be restored next time.
   if(is_local) {
     # !!! n_steps and sdr_assay fixed for now
     save_parameters(smoothing_fraction,
                     g_values,
                     extinction_threshold,
                     peak_split_depth,
                     sdr_coef,
                     nearest_nadir_distance)
   }
                 
                 
  param_list <- restore_param_list()
    
  # Returns NULL if the import fails. Returns a list ExperimentName and HormoneDF
  # if the import succeeds             
  input_folder_name <- v$data_folder_path
                 
  # Convert the files in the folder into a pulsar-formatted df
  input_data <- read_pulsar_input_folder(input_folder_name)
                 
  if (is.null(input_data)){
    showNotification("Error reading input folder", type="error")
  } else {
    
     not_id = showNotification("Processing folder...", duration = NULL)
                   
     v$input_data <- input_data
     v$experiment_name <- v$input_data$ExperimentName
     v$hormone_df <- v$input_data$HormoneDF
     v$all_results = NULL
     
     # !! Fixing the y-axes to the maximum value in the data columns (from column 3)
     if (input$chk_fix_y_axis){
       max_raw_conc <- max(v$hormone_df[, c(-1, -2)])
       } else {
         max_raw_conc = NULL
       }
                   
                   
    # Call the main pulsar computation routine, passing required args and storing result
    all_results = analyse_pulsar_data(v$hormone_df,
                                      smoothing_fraction,
                                      g_values,
                                      extinction_threshold,
                                      peak_split_depth,
                                      sdr_coef,
                                      nearest_nadir_distance,
                                      y_axis_max = max_raw_conc)
    if (!is.null(all_results)){
      v$all_results <- all_results
      
      
      clear_and_draw_plots(input, output, session, v$all_results)            
      
      # First tab is displayed by default. Load up here.
      # Other tabs load on tab click event
      result_index <- 1
      fill_tables(input, output, session, result_index, v$all_results)

      removeNotification(not_id)
      } else {
        # all_results came back NULL
        removeNotification(not_id)
        showNotification("Error when processing folder", type="error")
      }  # run_pulsar_batch returned null
                   
    } # end no error reading input folder
                 
  }) # end observer btn_batch
  

  #===================================================================
  #===================================================================

  # Save control shown depends on local vs. remote
  
#===================================================================
# On local save, write out current results organised into folders
observeEvent(input$btn_save,
     {
       not_id = showNotification("Saving output files...", duration = NULL)
                 
       # Current input control values were saved out on Run button
       # Read the saved file in here to pass to write method
       param_list <- restore_param_list()

       write_multi_pulsar(v$experiment_name, param_list, v$all_results, zip_files = FALSE)
                 
       removeNotification(not_id)
       showNotification("Finished saving...")
    }) #end observeEvent input$btn_save
  
  
#===================================================================
# On remote save, download current results as zip with dialogue 
  
output$btn_download <- downloadHandler(
  
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename <-   function() {
    paste(v$experiment_name, "_pulsar_" , get_current_time_string(), ".zip", sep="")
  } ,
  
  # This function should write data to a file given to it by
  # the argument 'file'.
  content = function(file)
  {
    not_id = showNotification("Saving output files...", duration = NULL)
    
    # Current input control values were saved out on Run button
    # Read the saved file in here to pass to write method
    param_list <- restore_param_list()
    
      
    # You're running from a web host...set the zip_file arg to true
      
    # They suggest writing to a temp directory on the server
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
      
    all_files <- write_multi_pulsar(v$experiment_name, param_list, v$all_results, zip_files = TRUE)
      
    removeNotification(not_id)
    showNotification("Finished saving...")
      
    zip_status <- zip::zip(file, all_files)
      
   #} # end not is_local

    # Delete files so they won't be downloaded again
    for (file in all_files)
    {
      unlink(file)
    }
  },
  contentType = "application/zip"
  
) # end downloadHandler
  
  
  
  
  
  
  # observeEvent(input$btn_save,
  #              {
  #                not_id = showNotification("Saving output files...", duration = NULL)
  #                
  #                # Current input control values were saved out on Run button
  #                # Read the saved file in here to pass to write method
  #                param_list <- restore_param_list()
  #                write_multi_pulsar(v$experiment_name, param_list, v$all_results)
  #                
  #                removeNotification(not_id)
  #                showNotification("Finished saving...")
  #              }) #end observeEvent input$btn_save
  
  
  
  
  
  
  
} # end server <- function()


#==========================
# Launch
#==========================

shinyApp(ui = ui, server = server)
