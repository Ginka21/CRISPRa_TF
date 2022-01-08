# 2021-12-27


# Load packages and source code -------------------------------------------

library("squash")
library("viridis")
library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "R_functions")
source(file.path(functions_dir, "01_calculating_scores.R"))
source(file.path(functions_dir, "02_labels_and_annotations.R"))



# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))




# Define functions --------------------------------------------------------

InvertTranspose <- function(input_mat) {
   t(apply(input_mat, 2, rev))
}

parula <- colorRampPalette(c("#352A87", "#0F5CDD", "#1481D6", "#06A4CA",
                             "#2EB7A4", "#87BF77", "#D1BB59", "#FEC832",
                             "#F9FB0E"
                             ))



MakeBreaks <- function(numeric_vec, num_breaks, trim = TRUE) {
   if (trim) {
      use_range <- quantile(numeric_vec, probs = c(0.02, 0.98), na.rm = TRUE)
   } else {
      use_range <- range(numeric_vec, na.rm = TRUE)
   }
   breaks_vec <- seq(from       = use_range[[1]],
                     to         = use_range[[2]],
                     length.out = num_breaks
                     )
   return(breaks_vec)
}


SigmoidalBreakPoints <- function(numeric_vec, use_for_log = 5, symm = TRUE) {

  numeric_vec <- as.numeric(numeric_vec)

  use_color_breaks <- exp(seq(log(use_for_log), log(100), length.out = 20)) - use_for_log
  use_color_breaks <- sort(unique(c(-use_color_breaks, use_color_breaks)))

  my_range <- range(numeric_vec, na.rm = TRUE)
  if (symm) {
    my_range <- c(-max(abs(my_range)), max(abs(my_range)))
  }
  my_range[1] <- my_range[1] - (diff(my_range) / 1000)
  my_range[2] <- my_range[2] + (diff(my_range) / 1000)

  use_color_breaks <- scales::rescale(use_color_breaks, to = my_range)
  use_color_breaks <- unique(round(use_color_breaks, digits = 3))
  return(use_color_breaks)
}



GetRepNumber <- function(column_name) {
   if (grepl("_rep", column_name)) {
      column_name_splits <- strsplit(column_name, "_", fixed = TRUE)[[1]]
      rep_string <- grep("^rep", column_name_splits, value = TRUE)
      rep_number <- as.integer(sub("rep", "", rep_string, fixed = TRUE))
      stopifnot(rep_number %in% 1:2)
   } else {
      rep_number <- NULL
   }
   return(rep_number)
}


BothRepColumns <- function(column_name) {
   rep_number <- GetRepNumber(column_name)
   if (rep_number == 1) {
      rep1_column <- column_name
      rep2_column <- sub("rep1", "rep2", column_name, fixed = TRUE)
   } else if (rep_number == 2) {
      rep1_column <- sub("rep2", "rep1", column_name, fixed = TRUE)
      rep2_column <- column_name
   }
   return(c(rep1_column, rep2_column))
}



BreaksForColumn <- function(input_df,
                            use_column,
                            num_empty_breaks          = 7,
                            num_other_breaks          = 50,
                            num_pos_breaks            = 10,
                            num_empty_to_other_breaks = 6,
                            num_other_to_pos_breaks   = 8,
                            use_custom_breaks         = NULL,
                            num_uniform_breaks        = NULL,
                            flatten_factor            = NULL,
                            weighting_for_controls    = NULL,
                            take_log2                 = FALSE
                            ) {


   if (!(is.null(use_custom_breaks))) {
      return(list(breaks = use_custom_breaks, type = "custom"))
   }

   mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
   are_empty <- input_df[, "Well_number_384"] %in% c(mat_384[, c(1, 24)])
   are_pos_ctrl <- input_df[, "Target_flag"] %in% "Pos. control"

   rep_number <- GetRepNumber(use_column)

   if (is.null(rep_number)) {
      all_vec <- input_df[, use_column]
   } else {
      rep_columns <- BothRepColumns(use_column)
      all_vec <- c(input_df[, rep_columns[[1]]], input_df[, rep_columns[[2]]])
      are_empty <- rep(are_empty, 2)
      are_pos_ctrl <- rep(are_pos_ctrl, 2)
   }

   if (IsPValue(use_column)) {
      all_vec <- -log10(all_vec)
   } else if (take_log2) {
      all_vec <- log2(all_vec)
   }

   use_diverging <- any(all_vec < 0)
   if (use_diverging) {
      most_extreme_value <- abs(range(all_vec))
      limits_vec <- c(-most_extreme_value, most_extreme_value)
   } else {
      limits_vec <- range(all_vec)
   }

   if (is.null(num_uniform_breaks) && use_diverging) {
      if (is.null(flatten_factor)) {
         if (grepl("SSMD", use_column, fixed = TRUE)) {
            flatten_factor <- 0.5
         } else {
            flatten_factor <- 1
         }
      }
      sigmoidal_breaks <- SigmoidalBreakPoints(limits_vec, use_for_log = flatten_factor)
      return(list(breaks = sigmoidal_breaks, type = "sigmoidal"))
   }

   if (is.null(weighting_for_controls)) {
      weighting_for_controls <- NeedsWeighting(use_column)
   }
   if (!(weighting_for_controls)) {
      num_uniform_breaks <- 100
   }

   if (!(is.null(num_uniform_breaks))) {
      uniform_breaks <- MakeBreaks(limits_vec, trim = FALSE, num_breaks = num_uniform_breaks)
      return(list(breaks = uniform_breaks, type = "uniform"))
   }


   if (any(all_vec < 0)) {
      most_extreme_value <- abs(range(all_vec))
      limits_vec <- c(-most_extreme_value, most_extreme_value)
      uniform_breaks <- SigmoidalBreakPoints(1:1000, use_for_log = 0.5)
      return(uniform_breaks)
   }


   empty_vec <- all_vec[are_empty]
   pos_vec   <- all_vec[are_pos_ctrl]
   other_vec <- all_vec[!(are_empty | are_pos_ctrl)]

   # stopifnot(max(empty_vec) < min(other_vec))
   # stopifnot(max(other_vec) < min(pos_vec))

   empty_breaks <- MakeBreaks(empty_vec, num_empty_breaks)
   empty_breaks <- c(min(all_vec), empty_breaks)

   other_breaks <- MakeBreaks(other_vec, num_other_breaks)

   pos_breaks <- MakeBreaks(pos_vec, num_pos_breaks)
   pos_breaks <- c(pos_breaks, max(all_vec))

   empty_to_other_breaks <- MakeBreaks(c(max(empty_breaks), min(other_breaks)),
                                       num_breaks = num_empty_to_other_breaks,
                                       trim = FALSE
                                       )

   other_to_pos_breaks   <- MakeBreaks(c(max(other_breaks), min(pos_breaks)),
                                       num_breaks = num_other_to_pos_breaks,
                                       trim = FALSE
                                       )

   all_breaks <- c(empty_breaks,
                   empty_to_other_breaks[2:(length(empty_to_other_breaks) - 1)],
                   other_breaks,
                   other_to_pos_breaks[2:(length(other_to_pos_breaks) - 1)],
                   pos_breaks
                   )
   all_breaks <- sort(unique(all_breaks))

   return(list(breaks = all_breaks, type = "weighted by controls"))
}


IsPValue <- function(column_name) {
   grepl("p_val", column_name, ignore.case = TRUE)
}

NeedsWeighting <- function(column_name) {
   grepl("_rep|Hit_strength", column_name)
}


PrettyRound <- function(numeric_vec) {
   old_scipen <- options("scipen" = 2)
   results_vec <- ifelse(abs(numeric_vec) > 100,
                         round(numeric_vec),
                         ifelse(abs(numeric_vec) > 10,
                                round(numeric_vec, digits = 1),
                                ifelse(abs(numeric_vec) > 1,
                                       round(numeric_vec, digits = 2),
                                       ifelse(abs(numeric_vec) > 0.01,
                                              signif(numeric_vec, digits = 2),
                                              signif(numeric_vec, digits = 1)
                                              )
                                       )
                                )
                         )
   results_vec <- as.character(results_vec)
   options(old_scipen)
   return(results_vec)
}


# lseq <- function(from, to, length.out) {
#   # logarithmic spaced sequence, from https://stackoverflow.com/a/29963530
#   exp(seq(log(from), log(to), length.out = length.out))
# }


HeatMap384 <- function(numeric_vec,
                       use_breaks,
                       ColorFunction  = NULL,
                       main_title     = "",
                       use_subtext    = "",
                       label_values   = FALSE,
                       use_minuslog10 = FALSE,
                       take_log2      = FALSE,
                       uniform_legend = FALSE
                       ) {

   assign("delete_numeric_vec", numeric_vec, envir = globalenv())
   assign("delete_use_breaks", use_breaks, envir = globalenv())

   if (is.null(ColorFunction)) {
      if (any(use_breaks < 0)) {
         ColorFunction <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
      } else {
         ColorFunction <- cividis
      }
   }

   stopifnot(length(numeric_vec) == 384)

   if (use_minuslog10) {
      final_vec <- -log10(numeric_vec)
   } else if (take_log2) {
      final_vec <- log2(numeric_vec)
   } else {
      final_vec <- numeric_vec
   }
   use_mat <- matrix(final_vec, nrow = 16, ncol = 24, byrow = TRUE)

   ## Create a matrix of colors for the heatmap
   use_mat <- matrix(final_vec, nrow = 16, ncol = 24, byrow = TRUE)
   use_cmap <- makecmap(use_mat, breaks = use_breaks, colFn = ColorFunction,
                        include.lowest = TRUE
                        )
   color_mat <- cmap(use_mat, map = use_cmap)


   ## Create a vector of colors for the legend
   if (uniform_legend) {
      num_steps <- 200
      legend_colors_seq <- seq(from       = min(use_breaks),
                               to         = max(use_breaks),
                               length.out = num_steps
                               )
      legend_colors_vec <- cmap(legend_colors_seq, map = use_cmap)
   } else {
      breaks_mat <- cbind(use_breaks[seq_len(length(use_breaks) - 1)],
                          use_breaks[seq_len(length(use_breaks) - 1) + 1L]
                          )
      legend_values_vec <- rowMeans(breaks_mat)
      legend_colors_vec <- use_cmap[["colors"]]
      num_steps <- length(legend_values_vec)
      stopifnot(length(legend_colors_vec) == num_steps)
   }

   ## Plot heatmap
   old_mai <- par(mai = c(1.3, 0.7, 1.1, 0.7))
   cimage(zcol = InvertTranspose(color_mat),
          xlab    = "",
          ylab    = "",
          xlabels = rep("", ncol(color_mat)),
          ylabels = rep("", nrow(color_mat)),
          tcl     = 0,
          mgp     = c(3, 0.3, 0),
          bty     = "n",
          axes    = FALSE
          )
   horizontal_lines <- seq(from = 0.5, to = nrow(use_mat) + 0.5, by = 1)
   vertical_lines   <- seq(from = 0.5, to = ncol(use_mat) + 0.5, by = 1)
   use_color <- "white"
   use_lwd <- 0.2
   segments(x0 = 0.5, x1 = ncol(use_mat) + 0.5, y0 = horizontal_lines,
            col = use_color, lwd = use_lwd
            )
   segments(x0 = vertical_lines, y0 = 0.5, y1 = nrow(use_mat) + 0.5,
            col = use_color, lwd = use_lwd
            )

   ## Plot heatmap legend

   step_width <- (par("usr")[[2]] - par("usr")[[1]]) / num_steps
   start_y <- par("usr")[[3]] - diff(grconvertY(c(0, 2), from = "lines", to = "user"))
   legend_height <- diff(grconvertX(c(0, 0.8), from = "lines", to = "user"))
   step_left_vec <- seq(from = par("usr")[[1]],
                        to   = par("usr")[[2]] - step_width,
                        by   = step_width
                        )
   rect(xleft    = step_left_vec,
        xright   = step_left_vec + step_width,
        ybottom  = start_y,
        ytop     = start_y + legend_height,
        col      = legend_colors_vec,
        border   = NA,
        xpd      = NA
        )
   if (uniform_legend) {
      assign("delete_legend_colors_sea", legend_colors_seq, envir = globalenv())
      legend_text_values <- pretty(legend_colors_seq, n = 10)
      assign("delete_legend_text_values_1", legend_colors_seq, envir = globalenv())
      use_tolerance <- diff(range(legend_colors_seq)) * 0.003
      are_within_limits <- (legend_text_values >= (min(legend_colors_seq) - use_tolerance)) &
                           (legend_text_values <= (max(legend_colors_seq) + use_tolerance))
      assign("delete_are_within_limits", are_within_limits, envir = globalenv())

      legend_text_values <- legend_text_values[are_within_limits]
      assign("delete_legend_text_values_1", legend_text_values, envir = globalenv())

      legend_text_vec <- format(legend_text_values, trim = TRUE)
      assign("delete_legend_text_vec", legend_text_vec, envir = globalenv())

      legend_text_range <- c(par("usr")[[1]] + (step_width / 2),
                             par("usr")[[2]] - (step_width / 2)
                             )
      assign("delete_legend_text_range", legend_text_range, envir = globalenv())
      legend_text_positions <- scales::rescale(legend_text_values,
                                               from = range(use_breaks),
                                               to   = legend_text_range
                                               )
   } else {
      use_indices <- seq(from = 1,
                         to = num_steps,
                         by = round((num_steps / 10) - 1)
                         )
      legend_text_values <- legend_values_vec[use_indices]
      legend_text_vec <- PrettyRound(legend_text_values)
      legend_text_positions <- (step_left_vec + (step_width / 2))[use_indices]
   }

   assign("delete_legend_text_positions", legend_text_positions, envir = globalenv())
   assign("delete_legend_text_value", legend_text_values, envir = globalenv())
   text(x      = legend_text_positions,
        y      = start_y - diff(grconvertY(c(0, 1.0), from = "lines", to = "user")),
        labels = legend_text_vec,
        adj    = c(0.5, 0),
        cex    = 0.8,
        xpd    = NA
        )
   border_grey <- "black"
   rect(xleft = par("usr")[[1]],
        xright = par("usr")[[2]],
        ybottom = start_y,
        ytop = start_y + legend_height,
        border = border_grey,
        xpd = NA
        )
   segments(x0 = legend_text_positions,
            y0 = start_y,
            y1 = start_y - diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
            col = border_grey,
            xpd = NA
            )
   ## Annotate plot
   text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 0.8), from = "lines", to = "user")),
        y      = seq_len(nrow(use_mat)),
        labels = rev(LETTERS[seq_len(nrow(use_mat))]),
        cex    = par("cex") * 0.8,
        col    = "black",
        xpd    = NA
        )
   text(x      = seq_len(ncol(use_mat)),
        y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
        labels = seq_len(ncol(use_mat)),
        cex    = par("cex") * 0.8,
        col    = "black",
        xpd    = NA
        )

   title(main_title, line = 2.9)
   mtext(use_subtext, side = 1, line = 4.2)


   ## Show values
   if (label_values) {
      x_positions <- rep(1:24, times = 16)
      y_positions <- rep(16:1, each = 24)
      use_dark_text <- colMeans(col2rgb(as.character(t(color_mat))) / 255) > 0.5
      text(x      = x_positions,
           y      = y_positions,
           labels = PrettyRound(numeric_vec),
           col    = ifelse(use_dark_text, "gray25", "gray75"),
           cex    = 0.4,
           font   = 2,
           xpd    = NA
           )
   }

   par(old_mai)
   return(invisible(NULL))
}



AveragedHeatmap <- function(input_df,
                            use_column,
                            both_replicates = TRUE,
                            main_title      = "",
                            use_subtext     = "",
                            ColorFunction   = NULL,
                            label_values    = FALSE,
                            take_log2       = FALSE,
                            uniform_legend  = NULL,
                            use_one_scale   = TRUE,
                            ...
                            ) {

   plates_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
   rep_number <- GetRepNumber(use_column)
   if (both_replicates && (!(is.null(rep_number)))) {
      rep_columns <- BothRepColumns(use_column)
      vec_list <- c(split(input_df[, rep_columns[[1]]], plates_vec),
                    split(input_df[, rep_columns[[2]]], plates_vec)
                    )
   } else {
      vec_list <- split(input_df[, use_column], plates_vec)
   }
   all_mat <- do.call(cbind, vec_list)
   mean_vec <- rowMeans(all_mat)

   if (!(use_one_scale) && both_replicates) {
      warning("The 'use_one_scale' parameter has no effect if 'both_replicates' is TRUE!")
   }

   if (use_one_scale || both_replicates) {
      breaks_list <- BreaksForColumn(input_df, use_column, take_log2 = take_log2, ...)
   } else {
      for (column_name in BothRepColumns(use_column)) {
         input_df[, column_name] <- input_df[, use_column]
      }
      breaks_list <- BreaksForColumn(input_df, use_column, take_log2 = take_log2, ...)
   }

   uniform_legend <- breaks_list[["type"]] == "uniform"

   HeatMap384(numeric_vec    = mean_vec,
              use_breaks     = breaks_list[["breaks"]],
              ColorFunction  = ColorFunction,
              main_title     = main_title,
              use_subtext    = use_subtext,
              label_values   = label_values,
              use_minuslog10 = IsPValue(use_column),
              take_log2      = take_log2,
              uniform_legend = uniform_legend
              )

   return(invisible(NULL))
}




HeatmapForPlate <- function(input_df,
                            plate_number,
                            use_column,
                            main_title     = NULL,
                            use_subtext    = "",
                            ColorFunction  = NULL,
                            show_z_prime   = NULL,
                            label_values   = FALSE,
                            take_log2      = FALSE,
                            uniform_legend = NULL,
                            use_one_scale  = TRUE,
                            ...
                            ) {

   has_replicates <- grepl("_rep", use_column, fixed = TRUE)

   if (is.numeric(plate_number)) {
      plate_number <- as.character(as.roman(plate_number))
   }
   are_this_plate <- input_df[, "Plate_number_384"] == plate_number

   if (is.null(main_title)) {
      main_title <- paste0("Plate ", plate_number)
      if (has_replicates) {
         rep_number <- GetRepNumber(use_column)
         main_title <- paste0(main_title, ", replicate ", rep_number)
      }
   }

   if (is.null(show_z_prime)) {
      show_z_prime <- has_replicates
   }

   if (show_z_prime) {
      z_prime <- Calculate_Z_Prime(input_df[are_this_plate, ], use_column)
      z_prime_text <- paste0("Z' = ", formatC(z_prime, digits = 2, format = "f"))
      if (is.null(use_subtext) || (nchar(use_subtext) == 0)) {
         use_subtext <- z_prime_text
      } else {
         use_subtext <- paste0(use_subtext, ", ", z_prime_text)
      }
   }

   if (use_one_scale) {
      breaks_list <- BreaksForColumn(input_df, use_column, take_log2 = take_log2, ...)
   } else {
      for (column_name in BothRepColumns(use_column)) {
         input_df[, column_name] <- input_df[, use_column]
      }
      breaks_list <- BreaksForColumn(input_df[are_this_plate, ], use_column,
                                     take_log2 = take_log2, ...
                                     )
   }
   if (is.null(uniform_legend)) {
      uniform_legend <- breaks_list[["type"]] == "uniform"
   }

   HeatMap384(numeric_vec    = input_df[are_this_plate, use_column],
              use_breaks     = breaks_list[["breaks"]],
              ColorFunction  = ColorFunction,
              main_title     = main_title,
              use_subtext    = use_subtext,
              label_values   = label_values,
              use_minuslog10 = IsPValue(use_column),
              take_log2      = take_log2,
              uniform_legend = uniform_legend
              )

   return(invisible(NULL))
}





# Draw example heatmaps ---------------------------------------------------

HeatmapForPlate(GBA_df, 10, "Raw_rep1", weighting_for_controls = FALSE,
                use_one_scale = FALSE
                )

HeatmapForPlate(GBA_df, 10, "Raw_rep1", weighting_for_controls = FALSE,
                use_one_scale = TRUE
                )


AveragedHeatmap(GBA_df, "Raw_rep1")
AveragedHeatmap(GBA_df, "Raw_rep1", both_replicates = FALSE)
AveragedHeatmap(GBA_df, "Raw_rep1", both_replicates = FALSE, use_one_scale = FALSE)
AveragedHeatmap(GBA_df, "Raw_rep2", both_replicates = FALSE)
AveragedHeatmap(GBA_df, "Raw_rep2", both_replicates = FALSE, use_one_scale = FALSE)

HeatmapForPlate(GBA_df, 1, "Raw_rep1")
HeatmapForPlate(GBA_df, 1, "Raw_rep1", use_one_scale = FALSE)
HeatmapForPlate(GBA_df, 1, "Raw_rep1", use_one_scale = TRUE, weighting_for_controls = FALSE)
HeatmapForPlate(GBA_df, 1, "Raw_rep1", use_one_scale = FALSE, weighting_for_controls = FALSE)



absolute_custom_breaks <- c(0, 250, seq(500, 1000, by = 50), 1500, seq(2000, 2500, by = 250))

HeatmapForPlate(GBA_df, 1, "Raw_rep1",
                use_custom_breaks = absolute_custom_breaks
                )
HeatmapForPlate(GBA_df, 1, "Raw_rep1",
                num_uniform_breaks = 50,
                ColorFunction = rocket
                )


HeatmapForPlate(GBA_df, 10, "Raw_rep2")


HeatmapForPlate(GBA_df, 2, "CellTiterGlo_raw")
HeatmapForPlate(GBA_df, 2, "CellTiterGlo_raw", uniform_legend = FALSE)


HeatmapForPlate(GBA_df, 2, "CellTiterGlo_foldNT")

HeatmapForPlate(GBA_df, 2, "CellTiterGlo_foldNT",
                weighting_for_controls = TRUE
                )


HeatmapForPlate(GBA_df, 2, "CellTiterGlo_foldNT",
                num_uniform_breaks = 100
                )

HeatmapForPlate(GBA_df, 12, "CellTiterGlo_raw",
                num_other_breaks = 100,
                )


HeatmapForPlate(GBA_df, 1, "SSMD_deltaNT")
HeatmapForPlate(GBA_df, 1, "SSMD_log2")

HeatmapForPlate(GBA_df, 1, "SSMD_deltaNT_Glo", uniform_legend = TRUE)



HeatmapForPlate(GBA_df, 6, "SSMD_deltaNT")



HeatmapForPlate(GBA_df, 10, "SSMD_log2")
HeatmapForPlate(GBA_df, 2, "p_value_log2", num_uniform_breaks = 100,
                use_subtext = expression(italic("p") * " value")
                )

HeatmapForPlate(GBA_df, 1, "Raw_rep1")
HeatmapForPlate(GBA_df, 2, "Raw_rep2")


HeatmapForPlate(GBA_df, 1, "DeltaNT_rep1")
HeatmapForPlate(GBA_df, 1, "DeltaNT_rep2")


HeatmapForPlate(GBA_df, 10, "Raw_rep1")
HeatmapForPlate(GBA_df, 10, "Raw_rep2")


HeatmapForPlate(GBA_df, 1, "SSMD_log2", label_values = TRUE)
HeatmapForPlate(GBA_df, 1, "Raw_rep1", label_values = TRUE)


HeatmapForPlate(GBA_df, 1, "FoldNT_rep1")
HeatmapForPlate(GBA_df, 1, "FoldNT_rep1", take_log2 = TRUE)


HeatmapForPlate(GBA_df, 5, "Hit_strength_deltaNT")



HeatmapForPlate(GBA_df, 10, "Raw_rep1", uniform_legend = TRUE)
HeatmapForPlate(GBA_df, 10, "Raw_rep1", uniform_legend = FALSE)

HeatmapForPlate(GBA_df, 10, "Raw_rep2")


HeatmapForPlate(GBA_df, 1, "SSMD_log2", label_values = TRUE)
HeatmapForPlate(GBA_df, 1, "Raw_rep1", label_values = TRUE)


HeatmapForPlate(GBA_df, 1, "FoldNT_rep1")
HeatmapForPlate(GBA_df, 1, "FoldNT_rep1", take_log2 = TRUE)

HeatmapForPlate(GBA_df, 1, "Raw_log2_rep1")
HeatmapForPlate(GBA_df, 1, "Raw_log2_rep1", weighting_for_controls = FALSE)



# Export heatmaps as PDF or PNG files -------------------------------------

heatmap_width <- 8
heatmap_height <- 6.5


for (label_cells in c(FALSE, TRUE)) {

   if (label_cells) {
      heatmaps_folder <- "Heatmaps (labelled)"
   } else {
      heatmaps_folder <- "Heatmaps"
   }

   for (i in seq_along(column_labels)) {

      current_column <- names(column_labels)[[i]]
      has_replicates <- !(is.null(GetRepNumber(current_column)))

      file_name <- paste0("Heatmaps - ", i, ") ", column_labels[[i]], ".pdf")
      file_name <- gsub("%", "percent", file_name, fixed = TRUE)
      column_subtext <- column_labels[[i]]

      pdf(file = file.path(output_dir, "Figures", heatmaps_folder, file_name),
          width = heatmap_width, height = heatmap_height
          )
      AveragedHeatmap(GBA_df, current_column,
                      main_title = "Mean of all plates and replicates (well effect)",
                      use_subtext = column_subtext,
                      label_values = label_cells
                      )
      for (plate_number in 1:12) {
         if (has_replicates) {
            rep_columns <- BothRepColumns(current_column)
            for (rep_column in rep_columns) {
               HeatmapForPlate(GBA_df, plate_number, rep_column,
                               use_subtext = column_subtext,
                               label_values = label_cells
                               )
            }
         } else {
            HeatmapForPlate(GBA_df, plate_number, current_column,
                            use_subtext = column_subtext,
                            label_values = label_cells
                            )
         }

      }
      dev.off()
   }

   for (i in seq_along(column_labels)) {

      current_column <- names(column_labels)[[i]]
      has_replicates <- !(is.null(GetRepNumber(current_column)))
      column_subtext <- column_labels[[i]]

      folder_name <- paste0("Heatmap PNGs - ", i, ") ", column_labels[[i]])
      folder_name <- gsub("%", "percent", folder_name, fixed = TRUE)
      folder_path <- file.path(output_dir, "Figures", heatmaps_folder, "PNGs", folder_name)
      dir.create(folder_path, showWarnings = FALSE)

      file_name <- paste0("Heatmap - ", column_labels[[i]],
                          " -- well effect"
                          )
      file_name <- gsub("%", "percent", file_name, fixed = TRUE)

      png(file = file.path(folder_path, file_name),
          width = heatmap_width, height = heatmap_height,
          units = "in", res = 600
          )
      AveragedHeatmap(GBA_df, current_column,
                      main_title = "Mean of all plates and replicates (well effect)",
                      use_subtext = column_subtext,
                      label_values = label_cells
                      )
      dev.off()
      for (plate_number in 1:12) {
         if (has_replicates) {
            rep_columns <- BothRepColumns(current_column)
            for (j in 1:2) {
               file_name <- paste0("Heatmap - ", column_labels[[i]],
                                   " - plate ", plate_number, " rep ", j
                                   )
               file_name <- gsub("%", "percent", file_name, fixed = TRUE)
               png(file = file.path(folder_path, file_name),
                   width = heatmap_width, height = heatmap_height,
                   units = "in", res = 600
                   )
               HeatmapForPlate(GBA_df, plate_number, rep_columns[[j]],
                               use_subtext = column_subtext,
                               label_values = label_cells
                               )
               dev.off()
            }
         } else {
            file_name <- paste("Heatmap - ", column_labels[[i]],
                               " - plate ", plate_number
                               )
            file_name <- gsub("%", "percent", file_name, fixed = TRUE)
            png(file = file.path(folder_path, file_name),
                width = heatmap_width, height = heatmap_height,
                units = "in", res = 600
                )
            HeatmapForPlate(GBA_df, plate_number, current_column,
                            use_subtext = column_subtext,
                            label_values = label_cells
                            )
            dev.off()
         }
      }
   }

}






