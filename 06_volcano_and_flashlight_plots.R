# 2121-12-27

# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/Screen"
functions_dir <- file.path(project_dir, "1_R_functions")
source(file.path(functions_dir, "03_plotting_helper_functions.R"))




# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))



# Define Functions --------------------------------------------------------

VolcanoFlashPlot <- function(input_df,
                             fc_column,
                             y_column,
                             show_only_genes = FALSE,
                             show_title = ""
                             ) {

   if (grepl("_rep", fc_column, fixed = TRUE)) {
      rep2_column <- sub("_rep1", "_rep2", fc_column, fixed = TRUE)
      fc_vec <- rowMeans(input_df[, c(fc_column, rep2_column)])
   } else {
      fc_vec <- input_df[, fc_column]
   }
   log_fc_vec <- log2(fc_vec)

   y_value_vec <- input_df[, y_column]
   if(grepl("SSMD_", y_column, fixed = TRUE)){
      y_label <- "SSMD"
   } else {
      y_value_vec <- -log10(y_value_vec)
      y_label <-"-log10 P value"
   }
   are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
   are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
   are_gene    <- !(is.na(input_df[, "Entrez_ID"]))


   if (show_only_genes) {
      are_valid <- are_gene
   } else {
      are_valid <- are_NT | are_posctrl | are_gene
   }

   use_margin <- c(4.25, 4, 3.5, 7.5)
   if (show_only_genes) {
      use_margin[[4]] <- 3.5
   }

   old_mar <- par(mar = use_margin)

   plot(1,
        xlim = range(log_fc_vec[are_valid]),
        ylim = range(c(0, y_value_vec[are_valid])),
        xlab = "log FC",
        ylab = y_label,
        las  = 1,
        mgp  = c(2.8, 0.7, 0),
        tcl  = -0.45,
        type = "n"
        )
   abline(h = 0, lty = "dotted", col = "grey70")
   abline(v = 0, lty = "dotted", col = "grey70")

   points(log_fc_vec[are_gene],
          y_value_vec[are_gene],
          pch = 16,
          col = adjustcolor("black", alpha.f = 0.3)
          )
   if (!(show_only_genes)) {

      pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
      NT_ctrl_color <- brewer.pal(5, "Blues")[[3]]

      points(log_fc_vec[are_posctrl],
             y_value_vec[are_posctrl],
             pch = 16,
             col = adjustcolor(pos_ctrl_color, alpha.f = 0.5)
             )

      points(log_fc_vec[are_NT],
             y_value_vec[are_NT],
             pch = 16,
             col = adjustcolor(NT_ctrl_color, alpha.f = 0.5)
             )

      controls_labels <- list(
         "NT"   = c("Non-targeting", "controls"),
         "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
         "Gene" = c("Genes in ", "CRISPRa", "library")
      )

      DrawSideLegend(labels_list = controls_labels,
                     use_colors = c(NT_ctrl_color, pos_ctrl_color, "black")
                     )
   }

   title(show_title, cex.main = 1.1)

   par(old_mar)

   return(invisible(NULL))
}



# Plot Data ---------------------------------------------------------------

VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value",
                 show_title = "Volcano Plot (p values from absolute deviations)"
                 )
VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value_log",
                 show_title = "Volcano Plot (p values from log fold change)"
                 )
VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo",
                 show_title = "Volcano Plot (standardized using CellTiterGlo)"
                 )
VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo_log",
                 show_title = "Volcano Plot (log FC, using CellTiterGlo)"
                 )


VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value", show_only_genes = TRUE)

VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_MM_paired", show_only_genes = TRUE,
                 show_title = "Dual-Flashlight Plot (GBA norm, SSMD from abs values)"
                 )

VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_log2", show_only_genes = FALSE,
                 show_title = "Dual-Flashlight Plot (GBA norm, SSMD log FC)"
                 )


# Export volcano plots as PDF and PNG -------------------------------------

base_width <- 5.5
base_height <- 5.1
use_dpi <- 600


for (only_genes in c(FALSE, TRUE)) {

   if (only_genes){
      file_name <- "Volcano plots - genes only.pdf"
      use_width <- base_width
   } else {
      file_name <- "Volcano plots - with controls.pdf"
      use_width <- base_width + 0.8
   }
   pdf(file = file.path(output_dir, "Figures", "Volcano_plots", file_name),
       width = use_width, height = base_height
       )

   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value",
               show_title = "Volcano Plot (p values from absolute deviations)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value_log",
               show_title = "Volcano Plot (p values from log fold change)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo",
               show_title = "Volcano Plot (standardized using CellTiterGlo)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo_log",
               show_title = "Volcano Plot (log FC, using CellTiterGlo)",
               show_only_genes = only_genes
               )
   dev.off()
}




for (only_genes in c(FALSE, TRUE)) {

   if (only_genes){
      file_postfix <- " - genes only.png"
      use_width <- base_width
   } else {
      file_postfix <- " - with controls.png"
      use_width <- base_width + 0.8
   }


   png(file = file.path(output_dir, "Figures", "Volcano_plots", "PNGs",
                        paste0("Volcano Plot - 1) p values from absolute deviations",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value",
               show_title = "Volcano Plot (p values from absolute deviations)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Volcano_plots", "PNGs",
                        paste0("Volcano Plot - 2) p values from log fold change",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "P_value_log",
               show_title = "Volcano Plot (p values from log fold change)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Volcano_plots", "PNGs",
                        paste0("Volcano Plot - 3) standardized using CellTiterGlo",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo",
               show_title = "Volcano Plot (standardized using CellTiterGlo)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Volcano_plots", "PNGs",
                        paste0("Volcano Plot - 4) log FC, using CellTiterGlo",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "P_value_Glo_log",
               show_title = "Volcano Plot (log FC, using CellTiterGlo)",
               show_only_genes = only_genes
               )
   dev.off()
}



# Export dual-flashlight plots as PDF and PNG -----------------------------


base_width <- 5.5
base_height <- 5.1
use_dpi <- 600


for (only_genes in c(FALSE, TRUE)) {

   if (only_genes){
      file_name <- "Dual-flashlight plots - genes only.pdf"
      use_width <- base_width
   } else {
      file_name <- "Dual-flashlight plots - with controls.pdf"
      use_width <- base_width + 0.8
   }
   pdf(file = file.path(output_dir, "Figures", "Dual_flashlight_plots", file_name),
       width = use_width, height = base_height
       )

   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_MM_paired",
               show_title = "Dual-flashlight Plot (Normalized data, SSMD MM)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_log2",
               show_title = "Dual-flashlight Plot  (Normalized data, SSMD MM log2)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "SSMD_MM_paired_Glo",
               show_title = "Dual-flashlight Plot (normalized using CellTiterGlo)",
               show_only_genes = only_genes
               )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "SSMD_log2_Glo",
               show_title = "Dual-flashlight Plot (normalized using CellTiterGlo, log2SSMD)",
               show_only_genes = only_genes
               )
   dev.off()
}




for (only_genes in c(FALSE, TRUE)) {

   if (only_genes){
      file_postfix <- " - genes only.png"
      use_width <- base_width
   } else {
      file_postfix <- " - with controls.png"
      use_width <- base_width + 0.8
   }


   png(file = file.path(output_dir, "Figures", "Dual_flashlight_plots", "PNGs",
                        paste0("Dual flashlight Plot - 1) SSMD MM paired",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_MM_paired",
               show_title = "Dual-flashlight Plot (Normalized data, SSMD MM)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Dual_flashlight_plots", "PNGs",
                        paste0("Dual-flashlight Plot - 2) Normalized data, SSMD MM log2",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_normalized", "SSMD_log2",
               show_title = "Dual-flashlight Plot (Normalized data, SSMD MM log2)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Dual_flashlight_plots", "PNGs",
                        paste0("Dual-flashlight Plot - 3) Glo normalized data, SSMD MM",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "SSMD_MM_paired_Glo",
               show_title = "Dual-flashlight Plot (Glo normalized data, SSMD MM)",
               show_only_genes = only_genes
               )
   dev.off()


   png(file = file.path(output_dir, "Figures", "Dual_flashlight_plots", "PNGs",
                        paste0("Dual-flashlight Plot - 4) log FC, using CellTiterGlo",
                               file_postfix
                               )
                        ),
       width = use_width, height = base_height, units = "in", res = use_dpi
       )
   VolcanoFlashPlot(GBA_df, "GBA_rep1_norm_Glo", "SSMD_log2_Glo",
               show_title = "Dual-flashlight Plot (Glo normalized, log2 SSMD)",
               show_only_genes = only_genes
               )
   dev.off()
}







