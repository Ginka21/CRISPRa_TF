# 2022-02-09


# Load packages and source code -------------------------------------------

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "07_Box_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir  <- file.path(project_dir, "3_R_objects", "2_GBA")
output_dir  <- file.path(project_dir, "4_output", "GBA")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Draw example plots ------------------------------------------------------

BeeBoxPlates(GBA_df, "Raw_rep1", plate_number = "II", split_NT = TRUE)

BeeBoxPlates(GBA_df, "Raw_rep1", show_subgroups = TRUE, plate_number = 2)
BeeBoxPlates(GBA_df, "Raw_rep1", show_subgroups = TRUE, show_96wp = TRUE)
BeeBoxPlates(GBA_df, "CellTiterGlo_raw", show_subgroups = TRUE, show_96wp = TRUE, plate_number = 22)
BeeBoxPlates(GBA_df, "CellTiterGlo_raw", show_subgroups = FALSE, show_96wp = TRUE, plate_number = 22)

BeeBoxPlates(GBA_df, "Raw_rep1", compare_group = "Gene")
BeeBoxPlates(GBA_df, "Raw_rep1", compare_group = "NT control")
BeeBoxPlates(GBA_df, "Raw_rep1", compare_group = "Pos. control")




# Export plots as PDF and PNG files ---------------------------------------

ExportAllBoxPlots(GBA_df, file.path(output_dir, "Figures", "Box plots"))








