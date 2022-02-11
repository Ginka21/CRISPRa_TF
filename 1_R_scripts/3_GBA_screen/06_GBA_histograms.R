#2022-01-04


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "04_Histograms.R"))



# Define folder path ------------------------------------------------------

r_data_dir  <- file.path(project_dir, "3_R_objects", "2_GBA")
output_dir  <- file.path(project_dir, "4_output", "GBA")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Plot histograms ---------------------------------------------------------

ThreeHistograms(GBA_df, "Raw_rep1")
ThreeHistograms(GBA_df, "FoldNT_rep1")
ThreeHistograms(GBA_df, "DeltaNT_rep1")
ThreeHistograms(GBA_df, "PercActivation_rep1")

ThreeHistograms(GBA_df, "Log2FC_rep1")


ThreeHistograms(GBA_df, "CellTiterGlo_raw")
ThreeHistograms(GBA_df, "CellTiterGlo_foldNT")
ThreeHistograms(GBA_df, "p_value_deltaNT")




# Export histograms as PDF and PNG files ----------------------------------

plot_width <- 7
plot_height <- 5


pdf(file = file.path(output_dir, "Figures", "Histograms", "Histograms.pdf"),
    width = plot_width, height = plot_height
    )
for (use_column in names(column_file_names)) {
  ThreeHistograms(GBA_df, use_column)
}
dev.off()



for (i in seq_along(column_file_names)) {
  use_column <- names(column_file_names)[[i]]
  file_name <- paste0("Histogram - ", i,  ") ", column_file_names[[i]], ".png")
  png(filename = file.path(output_dir, "Figures", "Histograms", "PNGs", file_name),
      width = plot_width, height = plot_height, units = "in", res = 600
      )
  ThreeHistograms(GBA_df, use_column)
  dev.off()
}


