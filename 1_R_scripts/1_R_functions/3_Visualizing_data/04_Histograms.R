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
    x_axis_label <- FormatPlotMath(long_column_labels[[use_column]])
  }

  ## Prepare data
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

  ## Determine x axis limits
  side_space_fraction <- 0.03
  data_range <- range(numeric_vec[are_valid])
  x_limits   <- DataAxisLimits(numeric_vec[are_valid], space_fraction = side_space_fraction)

  ## Compute the 3 histograms
  use_breaks <- seq(from = data_range[[1]], to = data_range[[2]], length.out = 70)
  gene_hist <- hist(numeric_vec[are_gene], breaks = use_breaks, plot = FALSE)
  NT_hist   <- hist(numeric_vec[are_NT],   breaks = use_breaks, plot = FALSE)
  pos_hist  <- hist(numeric_vec[are_pos],  breaks = use_breaks, plot = FALSE)

  ## Determine y axis limits
  count_max <- max(c(gene_hist[["counts"]], NT_hist[["counts"]], pos_hist[["counts"]]))
  y_limits <- c(count_max * -(side_space_fraction), count_max * (1 + side_space_fraction))

  ## Prepare graphical parameters
  fill_colors <- c(
    "gray55",                                            # genes
    colorRampPalette(brewer.pal(9, "Blues"))(100)[[80]], # NT control
    brewer.pal(5, "Reds")[[4]]                           # positive control
  )
  fill_colors <- vapply(fill_colors, function(x) adjustcolor(x, alpha.f = 0.8), "")
  border_colors <- c(
    "gray55",                    # genes
    brewer.pal(5, "Blues")[[5]], # NT control
    brewer.pal(5, "Reds")[[5]]   # positive control
  )

  use_mgp    <- c(2.7, 0.65, 0)
  use_tcl    <- -0.45

  ## Prepare the plot region
  old_mar    <- par(mar = c(5, 4.6, 3.8, 7.5))
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
  axis(1,          mgp = use_mgp, tcl = use_tcl)
  axis(2, las = 2, mgp = use_mgp, tcl = use_tcl)
  if (use_column == "CellTiterGlo_foldNT") {
    abline(v = 1, lty = "dashed", col = "gray40")
  }

  ## Add the 3 superimposed histograms
  if (grepl("p_?val", use_column, ignore.case = TRUE)) { # Prevent the positive controls from obscuring genes with low p values
    PlotHistResults(pos_hist,  fill_colors[[3]], border_colors[[3]])
    PlotHistResults(gene_hist, fill_colors[[1]], border_colors[[1]])
    PlotHistResults(NT_hist,   fill_colors[[2]], border_colors[[2]])
  } else {
    PlotHistResults(gene_hist, fill_colors[[1]], border_colors[[1]])
    PlotHistResults(NT_hist,   fill_colors[[2]], border_colors[[2]])
    PlotHistResults(pos_hist,  fill_colors[[3]], border_colors[[3]])
  }
  box(bty = "l")

  # Add legend
  controls_labels <- list(
    "Gene" = c("Genes in ", "CRISPRa", "library"),
    "NT"   = c("Non-targeting", "controls"),
    "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)"))
  )
  DrawSideLegend(labels_list    = controls_labels,
                 use_pch        = 22,
                 use_colors     = fill_colors,
                 border_colors  = border_colors,
                 use_point_size = 1.4,
                 lines_x_start  = 0.75
                 )

  return(invisible(NULL))
}





# Plot histograms ---------------------------------------------------------

ThreeHistograms(GBA_df, "Raw_rep1")
ThreeHistograms(GBA_df, "FoldNT_rep1")
ThreeHistograms(GBA_df, "DeltaNT_rep1")
ThreeHistograms(GBA_df, "PercActivation_rep1")

ThreeHistograms(GBA_df, "CellTiterGlo_raw")
ThreeHistograms(GBA_df, "CellTiterGlo_foldNT")
ThreeHistograms(GBA_df, "p_value_deltaNT")





# Export plots as PDF and PNG ---------------------------------------------

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




