% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulsar_utilities.R
\name{write_multi_pulsar}
\alias{write_multi_pulsar}
\title{write_multi_pulsar}
\usage{
write_multi_pulsar(
  experiment_name,
  param_list,
  pulsar_result_list,
  zip_files = FALSE
)
}
\arguments{
\item{experiment_name}{String value for experiment name. Typically one
uses the initialised value from reading and input file or folder.}

\item{param_list}{Named list containing algorithm parameters as produced by
restore_param_list}

\item{pulsar_result_list}{A list of pulsar results  (named fields $plot,
$peak_descriptors & $peak_summary) as produced by functions run_pulsar and
run_pulsar_batch.}

\item{zip_files}{Set to true to generate single zip containing all outfile png and csv.
Can be used for remote download when deployed as a shiny app.}
}
\description{
Generates output files from a complete pulsar analysis.
Accepts a list of pulsar results ($plot, $peak_descriptors & $peak_summary) as produced
by end-user-facing function run_pulsar and run_pulsar_batch.
Generates combined output files for all data series and individual images files
for each data series.
}
