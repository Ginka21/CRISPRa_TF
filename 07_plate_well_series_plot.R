# 2022-01-02


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_functions")
source(file.path(functions_dir, "02_labels_and_annotations.R"))
source(file.path(functions_dir, "03_plotting_helper_functions.R"))



# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))



# Define Functions --------------------------------------------------------

PlateWellPlot <- function(input_df,
                          use_column = "Raw_rep1",
                          show_title = "Plate-well series plot",
                          y_axis_label = NULL,
                          point_size = 0.6
                          ) {

   if (is.null(y_axis_label)) {
      if (("column_labels" %in% ls(envir = globalenv())) &&
          (use_column %in% names(column_labels))
          ) {
         y_axis_label <- column_labels[[use_column]]
      } else {
         y_axis_label <- use_column
      }
   }

   are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
   are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
   are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
   are_valid   <- are_NT | are_posctrl | are_gene

   plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
   has_replicates <- grepl("_rep", use_column, fixed = TRUE)

   if (has_replicates) {
      DoubleByPlate <- function(input_vec) {
         input_list <- split(input_vec, plate_numbers_vec)
         use_indices <- rep(seq_along(input_list), each = 2)
         doubled_list <- input_list[use_indices]
         doubled_vec <- unlist(doubled_list, use.names = FALSE)
         return(doubled_vec)
      }
      rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
      rep1_list <- split(input_df[, use_column], plate_numbers_vec)
      rep2_list <- split(input_df[, rep2_column], plate_numbers_vec)
      combined_list <- mapply(function(x, y) c(x, y), rep1_list, rep2_list, SIMPLIFY = FALSE)
      numeric_vec <- unlist(combined_list, use.names = FALSE)
      are_NT      <- DoubleByPlate(are_NT)
      are_posctrl <- DoubleByPlate(are_posctrl)
      are_gene    <- DoubleByPlate(are_gene)
      are_valid   <- DoubleByPlate(are_valid)
      plate_numbers_vec <- DoubleByPlate(plate_numbers_vec)
      replicate_numbers_vec <- rep(rep(1:2, each = 384), times = 12)
      plate_reps_vec <- paste0(plate_numbers_vec, "-rep", replicate_numbers_vec)
   } else {
      numeric_vec <- input_df[, use_column]
   }

   x_position_vec <- seq_along(numeric_vec) - 1L

   use_margin <- c(7, 4.6, 3.8, 7.5)
   old_mar <- par(mar = use_margin)
   use_mgp <- c(3, 0.65, 0)

   x_limits <- range(x_position_vec)
   x_limits <- x_limits + (diff(x_limits) * 0.04 * c(-1, 1))

   valid_vec <- numeric_vec[are_valid]

   # If the minimum and maximum are closer together than the distance
   # from the minimum to zero, then the y axis should not include zero
   far_from_zero <- diff(range(valid_vec)) < min(valid_vec)
   if (far_from_zero) {
      y_limits <- range(valid_vec)
   } else {
      y_limits <- range(c(0, valid_vec))
   }
   y_space <- diff(y_limits) * 0.04
   y_limits <- y_limits + (y_space * c(-1, 1))

   # If the minimum value is not negative or close to zero, the y range should stop at zero
   if (!(far_from_zero) && (min(valid_vec) >= y_space)) {
      y_limits[[1]] <- 0
   }

   plot(1,
        xlim = x_limits,
        ylim = y_limits,
        xaxs = "i",
        yaxs = "i",
        axes = FALSE,
        type = "n",
        ann = FALSE
        )
   box()
   title(show_title, cex.main = 1.1)
   axis(2, las = 1, mgp = use_mgp)
   axis(1, mgp = use_mgp, cex.axis = par("cex"))
   mtext(y_axis_label, side = 2, line = 3)
   label_x_pos <- par("usr")[[1]] + diff(grconvertX(c(0, 0.2), from = "lines", to = "user"))
   mtext("Well #:",
         side = 1,
         line = use_mgp[[2]],
         at   = label_x_pos,
         adj  = 1
         )

   plate_x_starts <- tapply(x_position_vec, plate_numbers_vec, function(x) x[[1]]) - 0.5
   plate_x_ends   <- tapply(x_position_vec, plate_numbers_vec, function(x) x[[length(x)]]) + 0.5

   rectangle_height <- diff(grconvertY(c(0, 0.9), from = "lines", to = "user"))
   start_y <- par("usr")[[3]] -
              diff(grconvertY(c(0, 3.2), from = "lines", to = "user"))
   text_colors <- c("gray25", "white")


   ## Draw the plate number indicator bar
   text(x      = label_x_pos,
        y      = start_y - (rectangle_height / 2),
        labels = "Plate #:", adj = c(1, 0.5),
        xpd    = NA
        )
   rect(xleft   = plate_x_starts,
        xright  = plate_x_ends,
        ybottom = start_y - rectangle_height,
        ytop    = start_y,
        col     = brewer.pal(9, "Blues")[c(3, 7)],
        xpd     = NA
        )
   text(x      = rowMeans(cbind(plate_x_starts, plate_x_ends)),
        y      = start_y - (rectangle_height / 2),
        labels = as.character(as.roman(unique(plate_numbers_vec))),
        cex    = par("cex") * 0.8,
        col    = text_colors,
        xpd    = NA
        )


   if (has_replicates) {

      rep_x_starts   <- tapply(x_position_vec, plate_reps_vec,    function(x) x[[1]]) - 0.5
      rep_x_ends   <- tapply(x_position_vec, plate_reps_vec,    function(x) x[[length(x)]]) + 0.5

      ## Draw the replicate number indicator bar
      rep_y_gap <- diff(grconvertY(c(0, 1.8), from = "lines", to = "user"))
      text(x      = label_x_pos,
           y      = start_y - rep_y_gap - (rectangle_height / 2),
           labels = "Replicate #:", adj = c(1, 0.5),
           xpd    = NA
           )
      rect(xleft   = rep_x_starts,
           xright  = rep_x_ends,
           ybottom = start_y - rep_y_gap - rectangle_height,
           ytop    = start_y - rep_y_gap,
           col     = brewer.pal(9, "Purples")[c(3, 7)],
           xpd     = NA
           )
      text(x      = rowMeans(cbind(rep_x_starts, rep_x_ends)),
           y      = start_y - rep_y_gap - (rectangle_height / 2),
           labels = rle(replicate_numbers_vec)[["values"]],
           cex    = par("cex") * 0.8,
           col    = text_colors,
           xpd    = NA
           )

   }

   points(x_position_vec[are_gene],
          numeric_vec[are_gene],
          pch = 16,
          col = adjustcolor("black", alpha.f = 0.3),
          cex = point_size
          )

   pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
   NT_ctrl_color <- brewer.pal(5, "Blues")[[3]]

   points(x_position_vec[are_posctrl],
          numeric_vec[are_posctrl],
          pch = 16,
          col = adjustcolor(pos_ctrl_color, alpha.f = 0.5),
          cex = point_size
          )

   points(x_position_vec[are_NT],
          numeric_vec[are_NT],
          pch = 16,
          col = adjustcolor(NT_ctrl_color, alpha.f = 0.5),
          cex = point_size
          )

   controls_labels <- list(
      "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
      "NT"   = c("Non-targeting", "controls"),
      "Gene" = c("Genes in ", "CRISPRa", "library")
   )

   DrawSideLegend(labels_list = controls_labels,
                  use_colors = c(pos_ctrl_color, NT_ctrl_color, "black")
                  )

   return(invisible(NULL))
}


# Draw Example Plots ------------------------------------------------------

PlateWellPlot(GBA_df)
PlateWellPlot(GBA_df, "FoldNT_rep1")
PlateWellPlot(GBA_df, "CellTiterGlo_raw")

PlateWellPlot(GBA_df, "DeltaNT_rep1")
PlateWellPlot(GBA_df, "Raw_log2_rep1")

PlateWellPlot(GBA_df, "Log2FC_rep1")

PlateWellPlot(GBA_df, "Hit_strength_deltaNT_Glo")



# Export Plots as PDF and PNG ---------------------------------------------

plot_width <- 7
plot_height <- 5.5

pdf(file = file.path(output_dir, "Figures", "Plate well series plots", "Plate well series plots.pdf"),
    width = plot_width, height = plot_height
    )
for (use_column in names(column_labels)) {
  PlateWellPlot(GBA_df, use_column, show_title = column_labels[[use_column]], y_axis_label = "")
}
dev.off()



for (i in seq_along(column_labels)) {
  use_column <- names(column_labels)[[i]]
  file_name <- paste0("Plate well series plot - ", i,  ") ",
                      sub("_rep1", "", use_column, fixed = TRUE),
                      ".png"
                      )
  png(file = file.path(output_dir, "Figures", "Plate well series plots", "PNGs", file_name),
      width = plot_width, height = plot_height, units = "in", res = 600
      )
  PlateWellPlot(GBA_df, use_column)
  dev.off()
}




