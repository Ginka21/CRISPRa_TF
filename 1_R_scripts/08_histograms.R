#2022-01-04


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "R_functions")
source(file.path(functions_dir, "01_calculating_scores.R"))
source(file.path(functions_dir, "02_labels_and_annotations.R"))
source(file.path(functions_dir, "03_plotting_helper_functions.R"))


# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")


# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))
plate_numbers_vec <- as.integer(as.roman(GBA_df[, "Plate_number_384"]))
split_df_list <- split(GBA_df, plate_numbers_vec)


# Define functions --------------------------------------------------------

PlotHistResults <- function(hist_results, fill_color, border_color) {
   half_width <- (hist_results[["mids"]][[2]] - hist_results[["mids"]][[1]]) / 2
   are_not_zero <- hist_results[["counts"]] != 0
   mids_vec <- hist_results[["mids"]][are_not_zero]
   rect(xleft   = mids_vec - half_width,
        xright  = mids_vec + half_width,
        ybottom = 0,
        ytop    = hist_results[["counts"]][are_not_zero],
        col     = fill_color,
        border  = border_color,
        lwd     = 0.5,
        xpd     = NA
        )
}

ThreeHistograms <- function(input_df, use_column, x_axis_label = NULL) {

   if (is.null(x_axis_label)) {
      x_axis_label <- column_labels[[use_column]]
   }

   are_gene  <- !(is.na(input_df[, "Entrez_ID"]))
   are_NT    <- input_df[, "Is_NT_ctrl"]
   are_pos   <- input_df[, "Is_pos_ctrl"]
   are_valid <- are_NT | are_pos | are_gene

   has_replicates <- grepl("_rep", use_column, fixed = TRUE)
   if (has_replicates) {
      rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
      numeric_vec <- rowMeans(input_df[, c(use_column, rep2_column)])
   } else {
      numeric_vec <- input_df[, use_column]
   }

   fill_colors <- c(
      "gray55",                                            # genes
      colorRampPalette(brewer.pal(9, "Blues"))(100)[[80]], # NT control
      brewer.pal(5, "Reds")[[4]]                           # positive_control
   )
   fill_colors <- vapply(fill_colors, function(x) adjustcolor(x, alpha.f = 0.8), "")
   border_colors <- c(
      "gray55",                                            # genes
      brewer.pal(5, "Blues")[[5]],                         # NT control
      brewer.pal(5, "Reds")[[5]]                           # positive_control
   )

   use_margin <- c(5, 4.6, 3.8, 7.5)
   old_mar    <- par(mar = use_margin)
   use_mgp    <- c(3, 0.65, 0)
   data_range <- range(numeric_vec[are_valid])
   x_limits   <- DataAxisLimits(numeric_vec[are_valid])

   use_breaks <- seq(from = data_range[[1]],
                     to = data_range[[2]],
                     length.out = 70
                     )

   # Plot First distribution
   hist_results <- hist(numeric_vec[are_gene],
                        breaks = use_breaks,
                        plot = FALSE
                        )
   count_max <- max(hist_results[["counts"]])
   y_limits <- c(count_max * -0.03, count_max * 1.03)

   assign("delete_hist_results", hist_results, envir = globalenv())

   use_mgp <- c(2.7, 0.65, 0)

   plot(1,
        xlim   = x_limits,
        ylim   = y_limits,
        xaxs   = "i",
        yaxs   = "i",
        axes   = FALSE,
        xlab   = x_axis_label,
        ylab   = "Count",
        mgp    = use_mgp,
        type   = "n"
        )

   use_tcl <- -0.45
   axis(1, mgp = use_mgp, tcl = use_tcl)
   axis(2, las = 2, mgp = use_mgp, tcl = use_tcl)

   PlotHistResults(hist_results, fill_colors[[1]], border_colors[[1]])

   hist_results <- hist(numeric_vec[are_NT], breaks = use_breaks, plot = FALSE)
   PlotHistResults(hist_results, fill_colors[[2]], border_colors[[2]])

   hist_results <- hist(numeric_vec[are_pos], breaks = use_breaks, plot = FALSE)
   PlotHistResults(hist_results, fill_colors[[3]], border_colors[[3]])

   box(bty = "l")

   # Add legend
   controls_labels <- list(
      "Gene" = c("Genes in ", "CRISPRa", "library"),
      "NT"   = c("Non-targeting", "controls"),
      "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)"))
   )

   DrawSideLegend(labels_list = controls_labels,
                  use_pch = 22,
                  use_colors = fill_colors,
                  border_colors = border_colors,
                  use_point_size = 1.4,
                  lines_x_start = 0.75
                  )

   return(invisible(NULL))
}





ThreeHistograms(GBA_df, "Raw_rep1")
ThreeHistograms(GBA_df, "FoldNT_rep1")
ThreeHistograms(GBA_df, "DeltaNT_rep1")
ThreeHistograms(GBA_df, "PercActivation_rep1")

ThreeHistograms(GBA_df, "CellTiterGlo_raw")
ThreeHistograms(GBA_df, "p_value_deltaNT")






# Export plots as PDF and PNG ---------------------------------------------

plot_width <- 5.5
plot_height <- 4

pdf(file = file.path(output_dir, "Figures", "Histograms", "Histograms.pdf"),
    width = plot_width, height = plot_height
    )
for (use_column in names(column_labels)) {
  ThreeHistograms(GBA_df, use_column)
}
dev.off()



for (i in seq_along(column_labels)) {
  use_column <- names(column_labels)[[i]]
  file_name <- paste0("Histogram - ", i,  ") ",
                      sub("_rep1", "", use_column, fixed = TRUE),
                      ".png"
                      )
  png(file = file.path(output_dir, "Figures", "Histograms", "PNGs", file_name),
      width = plot_width, height = plot_height, units = "in", res = 600
      )
  ThreeHistograms(GBA_df, use_column)
  dev.off()
}


