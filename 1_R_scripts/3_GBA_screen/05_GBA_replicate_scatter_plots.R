# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "03_Replicate_scatter_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "2_GBA")
output_dir <- file.path(project_dir, "4_output", "GBA")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Examine the correlation between replicates ------------------------------

ReplicateScatter(GBA_df, "Raw_rep1")
ReplicateScatter(GBA_df, "PercActivation_log2_Glo_rep1")

ReplicateScatter(GBA_df, "Log2FC_rep1", same_scale = FALSE)




# Export plots as PDF and PNG files ---------------------------------------

ExportAllReplicateScatterPlots(GBA_df)



