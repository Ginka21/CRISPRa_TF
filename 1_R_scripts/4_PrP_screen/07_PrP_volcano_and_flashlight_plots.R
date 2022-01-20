# 2022-01-20


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "05_Volcano_and_flashlight_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir,"4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))




# Define labels -----------------------------------------------------------

controls_labels[["Pos"]] <- c("Positive", "controls", expression("(" * italic("PRNP") * " gene)"))




# Plot data ---------------------------------------------------------------

VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_deltaNT",
                 show_title = "Volcano plot (p values from untransformed data)"
                 )

VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_act",
                 show_title = "Volcano plot (p values from % activation)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = "Volcano plot (p values from log2-transformed data)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_act_log2",
                 show_title = "Volcano plot (p values from % activation, log2 data)"
                 )

VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_deltaNT_Glo",
                 show_title = "Volcano plot (untransformed, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_act_Glo",
                 show_title = "Volcano plot (% activation, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = "Volcano plot (log2, CellTitreGlo-normalized)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_act_log2_Glo",
                 show_title = "Volcano plot (% activation, log2, CellTitreGlo-norm.)"
                 )


VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )

VolcanoFlashPlot(PrP_df, "PercActivation_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )




# Export individual plots as PNG files ------------------------------------

base_width <- 5.5
base_height <- 5.1

selected_volcanoes_dir <- file.path(output_dir, "Figures", "Volcano plots", "Selected plots")


png(file.path(selected_volcanoes_dir, "1) Volcano plot - cutoffs shown.png"),
    width = base_width + 0.8, height = base_height, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = Embolden(FormatPlotMath("Volcano plot (p values, log2FC)")),
                 label_points = FALSE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), indicate_p_values = 0.05
                 )
dev.off()


pdf(file.path(selected_volcanoes_dir, "2) Volcano plot - genes shown.pdf"),
    width = base_width + 0.8, height = base_height#, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = Embolden(FormatPlotMath("Volcano plot (p values, log2FC)")),
                 label_points = TRUE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), tiny_labels = TRUE
                 )
dev.off()


png(file.path(selected_volcanoes_dir, "3) Volcano plot - Glo-normalized - cutoffs shown.png"),
    width = base_width + 0.8, height = base_height, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = Embolden(FormatPlotMath("Volcano plot (normalized to CellTiter-Glo)")),
                 label_points = FALSE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2)
                 )
dev.off()


stop()

# Export plots as PDF and PNG files ---------------------------------------

ExportAllVolcanoAndFlashlightPlots(PrP_df)





