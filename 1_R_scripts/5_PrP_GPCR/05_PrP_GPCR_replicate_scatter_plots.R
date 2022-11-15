# 2022-09-19


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "03_Replicate_scatter_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP", "GPCRa")
output_dir <- file.path(project_dir, "4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Modify labels (for PrPc screen) -----------------------------------------

are_Glo <- grepl("Glo", names(column_file_names), fixed = TRUE)
column_file_names <- column_file_names[!(are_Glo)]
AdjustLabels()



# Examine the correlation between replicates ------------------------------

ReplicateScatter(PrP_df, "Raw_rep1")
ReplicateScatter(PrP_df, "PercActivation_log2_rep1")
ReplicateScatter(PrP_df, "Log2FC_rep1", same_scale = FALSE)





# Export plots as PDF and PNG files ---------------------------------------

ExportAllReplicateScatterPlots(PrP_df,
                               file.path(output_dir, "Figures", "GPCRa",
                                         "Replicate scatter plots"
                                         )
                               )



