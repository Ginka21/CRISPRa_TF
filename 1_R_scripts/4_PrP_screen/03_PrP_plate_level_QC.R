# 2022-01-18


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "2_Analyzing_data",    "01_calculating_scores.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "01_Plate_level_QC.R"))



# Define folder path ------------------------------------------------------

r_data_dir  <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir  <- file.path(project_dir,"4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))




# Calculate plate-wise quality metrics ------------------------------------

PlotZPrimes(PrP_df)
PlotSSMDControls(PrP_df)

plot_width <- 5.5
plot_height <- 3.8

pdf(file = file.path(output_dir, "Figures", "Quality metrics", "Quality metrics.pdf"),
    width = plot_width, height = plot_height
    )
PlotZPrimes(PrP_df)
PlotSSMDControls(PrP_df)
dev.off()


png(filename = file.path(output_dir, "Figures", "Quality metrics", "Z_prime.png"),
    width = plot_width, height = plot_height, units = "in", res = 600
    )
PlotZPrimes(PrP_df)
dev.off()


png(filename = file.path(output_dir, "Figures", "Quality metrics", "SSMD.png"),
    width = plot_width, height = plot_height, units = "in", res = 600
    )
PlotSSMDControls(PrP_df)
dev.off()


