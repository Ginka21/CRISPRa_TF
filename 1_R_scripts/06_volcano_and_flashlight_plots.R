# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "R_functions")
source(file.path(functions_dir, "03_plotting_helper_functions.R"))



# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))



# Define Functions --------------------------------------------------------

VolcanoFlashPlot <- function(input_df,
                             log_fc_column,
                             y_column,
                             show_only_genes = FALSE,
                             show_title = "",
                             point_size = 0.6
                             ) {

   if (grepl("_rep", log_fc_column, fixed = TRUE)) {
      rep2_column <- sub("_rep1", "_rep2", log_fc_column, fixed = TRUE)
      log_fc_vec <- rowMeans(input_df[, c(log_fc_column, rep2_column)])
   } else {
      log_fc_vec <- input_df[, log_fc_column]
   }
   y_value_vec <- input_df[, y_column]

   if (grepl("PercActivation", log_fc_column, fixed = TRUE)) {
      x_label <- "% activation"
      log_fc_vec <- log_fc_vec * 100
   } else {
      x_label <- expression("log"["2"] ~ "fold change")
   }
   if (grepl("SSMD", y_column, fixed = TRUE)) {
      y_label <- "SSMD"
   } else {
      y_value_vec <- -log10(y_value_vec)
      y_label <- expression("" - "log"["10"] ~ plain("p") ~ "value")
   }
   are_NT      <- input_df[, "Is_NT_ctrl"]
   are_posctrl <- input_df[, "Is_pos_ctrl"]
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
        xlab = x_label,
        ylab = y_label,
        las  = 1,
        mgp  = c(2.8, 0.7, 0),
        tcl  = -0.45,
        type = "n",
        bty  = "n"
        )
   abline(h = 0, lty = "dotted", col = "grey70")
   abline(v = 0, lty = "dotted", col = "grey70")
   box()

   points(log_fc_vec[are_gene],
          y_value_vec[are_gene],
          pch = 16,
          col = adjustcolor("black", alpha.f = 0.3),
          cex = point_size
          )
   if (!(show_only_genes)) {

      pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
      NT_ctrl_color <- brewer.pal(5, "Blues")[[3]]

      points(log_fc_vec[are_posctrl],
             y_value_vec[are_posctrl],
             pch = 16,
             col = adjustcolor(pos_ctrl_color, alpha.f = 0.5),
             cex = point_size
             )

      points(log_fc_vec[are_NT],
             y_value_vec[are_NT],
             pch = 16,
             col = adjustcolor(NT_ctrl_color, alpha.f = 0.5),
             cex = point_size
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



# Define pairs of metrics -------------------------------------------------

pairs_list <- list(

   "deltaNT" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_deltaNT",
      "title" = "Volcano plot (p values from untransformed data)"
   ),
   "activation" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_act",
      "title" = "Volcano plot (p values from % activation)"
   ),
   "log2" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_log2",
      "title" = "Volcano plot (p values from log2-transformed data)"
   ),
   "act_log2" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_act_log2",
      "title" = "Volcano plot (p values from % activation, log2 data)"
   ),


   "deltaNT_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_deltaNT_Glo",
      "title" = "Volcano plot (untransformed data, CellTitreGlo-norm.)"
   ),
   "activation_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_act_Glo",
      "title" = "Volcano plot (% activation, CellTitreGlo-norm.)"
   ),
   "log2_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_log2_Glo",
      "title" = "Volcano plot (log2, CellTitreGlo-normalized)"
   ),
   "act_log2_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_act_log2_Glo",
      "title" = "Volcano plot (% activation, log2, CellTitreGlo-norm.)"
   )


)


# Plot Data ---------------------------------------------------------------

VolcanoFlashPlot(GBA_df, "Log2FC_rep1", "p_value_deltaNT",
                 show_title = "Volcano plot (p values from untransformed data)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_rep1", "p_value_act",
                 show_title = "Volcano plot (p values from % activation)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_rep1", "p_value_log2",
                 show_title = "Volcano plot (p values from log2-transformed data)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_rep1", "p_value_act_log2",
                 show_title = "Volcano plot (p values from % activation, log2 data)"
                 )

VolcanoFlashPlot(GBA_df, "Log2FC_Glo_rep1", "p_value_deltaNT_Glo",
                 show_title = "Volcano plot (untransformed, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_Glo_rep1", "p_value_act_Glo",
                 show_title = "Volcano plot (% activation, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = "Volcano plot (log2, CellTitreGlo-normalized)"
                 )
VolcanoFlashPlot(GBA_df, "Log2FC_Glo_rep1", "p_value_act_log2_Glo",
                 show_title = "Volcano plot (% activation, log2, CellTitreGlo-norm.)"
                 )


VolcanoFlashPlot(GBA_df, "Log2FC_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )

VolcanoFlashPlot(GBA_df, "PercActivation_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )








# Export volcano plots as PDF and PNG -------------------------------------

base_width <- 5.5
base_height <- 5.1

plot_types <- c("Volcano", "Dual flashlight (logFC)", "Dual flashlight (% activation)")

for (use_device in c("none", "pdf", "png")) {

   for (plot_type in plot_types) {

      if (plot_type == "Volcano") {
         folder_name <- "Volcano_plots"
         PNG_prefix <- "Volcano plot"
         x_sub <- "Log2FC"
         y_sub <- "p_value"
         title_prefix <- "Volcano plot"
      } else if (plot_type == "Dual flashlight (logFC)") {
         folder_name <- "Dual_flashlight_plots_logFC"
         PNG_prefix <- "Dual flashlight (logFC)"
         x_sub <- "Log2FC"
         y_sub <- "SSMD"
         title_prefix <- "Dual flashlight plot"
      } else if (plot_type == "Dual flashlight (% activation)") {
         folder_name <- "Dual_flashlight_plots_percent_activation"
         PNG_prefix <- "Dual flashlight (% activation)"
         x_sub <- "PercActivation"
         y_sub <- "SSMD"
         title_prefix <- "Dual flashlight plot"
      }

      for (only_genes in c(FALSE, TRUE)) {

         file_name <- folder_name

         if (only_genes) {
            PDF_name <- paste0(file_name, " - genes only.pdf")
            use_width <- base_width
         } else {
            PDF_name <- paste0(file_name, " - with controls.pdf")
            use_width <- base_width + 0.8
         }
         if (use_device == "pdf") {
            pdf(file = file.path(output_dir, "Figures", folder_name, PDF_name),
                width = use_width, height = base_height
                )
         }
         for (i in seq_along(pairs_list)) {
            plot_title <- sub("Volcano plot", title_prefix, pairs_list[[i]][["title"]])
            if (plot_type != "Volcano") {
               plot_title <- sub("p values", "SSMD", plot_title, fixed = TRUE)
            }
            if (use_device == "png") {
               PNG_name <- paste0(i, ") ", gsub("%", "percent", plot_title, fixed = TRUE))
               if (only_genes) {
                  PNG_name <- paste0(PNG_name, " - genes only.png")
               } else {
                  PNG_name <- paste0(PNG_name, " - with controls.png")
               }
               png(file = file.path(output_dir, "Figures", folder_name, "PNGs", PNG_name),
                   width = use_width, height = base_height, units = "in", res = 600
                   )
            }
            VolcanoFlashPlot(GBA_df,
                             sub("Log2FC", x_sub, pairs_list[[i]][["x_var"]], fixed = TRUE),
                             sub("p_value", y_sub, pairs_list[[i]][["y_var"]], fixed = TRUE),
                             show_title = plot_title,
                             show_only_genes = only_genes
                             )
            if (use_device == "png") {
               dev.off()
            }
         }

         if (use_device == "pdf") {
            dev.off()
         }
      }
   }

}



