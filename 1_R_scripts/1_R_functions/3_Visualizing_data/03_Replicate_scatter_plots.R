# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

ScatterPlot <- function(x_vec,
                        y_vec,
                        use_limits      = NULL,
                        point_size      = 0.6,
                        label_y_axis    = TRUE,
                        use_color       = NULL,
                        top_label       = NULL,
                        draw_regression = TRUE
                        ) {

  stopifnot(length(x_vec) == length(y_vec))

  ## Determine axis limits
  if (is.null(use_limits)) {
    use_limits = range(c(x_vec, y_vec))
  }
  use_limits <- use_limits + (diff(use_limits) * 0.04 * c(-1, 1))

  if (draw_regression) {
    ## Perform a linear regression (and compute a 95% confidence interval)
    model_df <- data.frame("x_var" = x_vec, "y_var" = y_vec)
    lm_model <- lm(y_var ~ x_var, data = model_df)
    lm_summary <- summary(lm_model)
    new_seq <- seq(use_limits[[1]], use_limits[[2]], length.out = 200)
    new_df <- data.frame("x_var" = new_seq)
    conf_int_mat <- predict(lm_model,
                            newdata = new_df,
                            interval = "confidence",
                            level = 0.95
                            )
    corr_text <- bquote(italic("R") * ""^2  ~ "=" ~
                        .(format(round(lm_summary[["r.squared"]], digits = 2), nsmall = 2))
                        )
  } else {
    pearsons_r <- cor.test(x_vec, y_vec)[["estimate"]][[1]]
    corr_text <- bquote(italic("r")  ~ "=" ~
                        .(format(round(pearsons_r, digits = 2), nsmall = 2))
                        )
  }

  ## Define graphical parameters
  use_mgp <- c(2.8, 0.55, 0)
  use_tcl <- -0.35
  if (is.null(use_color)) {
    use_color <- "#000000"
  }

  ## Set up the plot region
  plot(1,
       xlim = use_limits,
       ylim = use_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       type = "n"
       )

  ## Draw and label axes
  axis(1, mgp = use_mgp, tcl = use_tcl, gap.axis = 0.5)
  mtext("Replicate 1", side = 1, line = 2.2)
  axis(2, las = 1, mgp = use_mgp, tcl = use_tcl)
  if (label_y_axis) {
    mtext("Replicate 2", side = 2, line = use_mgp[[1]])
  }

  ## Draw the plot title
  if (!(is.null(top_label))) {
    mtext(Embolden(VerticalAdjust(top_label)),
          line = 1.6, cex = par("cex"), font = 2
          )
  }
  mtext(VerticalAdjust(as.expression(corr_text)),
        line = 0.05, cex = par("cex"), font = 2
        )

  ## Draw indicator lines
  abline(a = 0, b = 1, col = "grey80", lty = "dashed")
  abline(h = 0, col = "gray90")
  abline(v = 0, col = "gray90")

  ## Draw the regression line and 95% CI region
  if (draw_regression) {
    polygon(c(new_df[, 1], rev(new_df[, 1])),
            c(conf_int_mat[, 2], rev(conf_int_mat[, 3])),
            col = Palify(use_color, fraction_pale = 0.8), border = NA
            )
    lines(new_df[, 1], conf_int_mat[, 1], col = use_color, lwd = 1.5)
  }
  box()

  ## Draw the points of the scatter plot
  points(x_vec,
         y_vec,
         pch = 16,
         col = adjustcolor(use_color, alpha.f = 0.5),
         cex = point_size * par("cex")
         )

  return(invisible(NULL))
}


ReplicateScatter <- function(input_df,
                             rep1_column,
                             rep2_column = NULL,
                             show_title = "Replicate scatter plot",
                             same_scale = TRUE,
                             ...
                             ) {

  if (is.null(rep2_column)) {
    rep2_column <- sub("rep1", "rep2", rep1_column, fixed = TRUE)
  }

  are_gene <- !(is.na(input_df[, "Entrez_ID"]))
  corr_gene <- cor.test(input_df[are_gene, rep1_column],
                        input_df[are_gene, rep2_column]
                        )[["estimate"]][[1]]
  are_NT      <- input_df[, "Is_NT_ctrl"]
  are_posctrl <- input_df[, "Is_pos_ctrl"]
  are_valid <- are_NT | are_posctrl | are_gene

  if (same_scale) {
    axis_limits <- range(input_df[are_valid, c(rep1_column, rep2_column)])
  } else {
    axis_limits <- NULL
  }

  ## Set up the multi-plot layout
  layout_mat <- rbind(c(1, 2, 2, 2, 2, 2, 3),
                      5:11,
                      rep(4, 7)
                      )
  use_heights <- c(0.35, 0.45, 0.2)
  use_widths <- c(0.11, 0.21, 0.1, 0.21, 0.1, 0.21, 0.06)
  layout(layout_mat,
         heights = use_heights,
         widths = use_widths
         )
  old_par <- par(mar = rep(0, 4), cex = par("cex") / 0.66)

  for (i in 1:2) {
    MakeEmptyPlot()
  }
  text(x = 0.5, y = 0.7, labels = show_title, cex = par("cex") * 1.1)
  for (i in 1:3) {
    MakeEmptyPlot()
  }

  ## Draw the 3 scatter plots
  pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
  NT_ctrl_color <- brewer.pal(5, "Blues")[[4]]

  ScatterPlot(input_df[are_gene, rep1_column],
              input_df[are_gene, rep2_column],
              top_label = "Transcription factors",
              use_limits = axis_limits,
              ...
              )
  MakeEmptyPlot()
  ScatterPlot(input_df[are_NT, rep1_column],
              input_df[are_NT, rep2_column],
              top_label = "NT controls",
              use_color = NT_ctrl_color,
              use_limits = axis_limits,
              label_y_axis = FALSE,
              ...
              )
  MakeEmptyPlot()
  ScatterPlot(input_df[are_posctrl, rep1_column],
              input_df[are_posctrl, rep2_column],
              top_label = "Positive controls",
              use_color = pos_ctrl_color,
              use_limits = axis_limits,
              label_y_axis = FALSE,
              ...
              )
  MakeEmptyPlot()

  par(old_par)
  layout(1)

  return(invisible(NULL))
}



PlotPlateQualities <- function(rep1_vec,
                               rep2_vec,
                               y_limits_include = NULL,
                               y_axis_label = "",
                               quality_ranges = list(c(0.5, 1),
                                                     c(0, 0.5),
                                                     c(-Inf, 0)
                                                     )
                               ) {

  stopifnot(length(rep1_vec) == length(rep2_vec))
  data_vec <- c(rep1_vec, rep2_vec)

  ## Prepare x axis positions
  x_mids <- seq_along(rep1_vec)
  x_space <- 0.5
  x_positions <- c(x_mids - (x_space / 2), x_mids + (x_space / 2))
  x_space <- 0.5 + length(x_mids) * 0.015
  x_limits <- c(1 - x_space, length(x_mids) + x_space)

  ## Prepare y axis positions
  y_limits <- range(c(y_limits_include, data_vec))
  y_space <- diff(y_limits) * 0.05
  if (y_limits[[1]] > (min(data_vec) - y_space)) {
    y_limits[[1]] <- y_limits[[1]] - y_space
  }
  if (y_limits[[2]] < (max(data_vec) + y_space)) {
    y_limits[[2]] <- y_limits[[2]] + y_space
  }

  ## Set up the plot region
  old_mai <- par("mai" = c(0.7, 0.82, 0.5, 0.42))
  plot(1,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  axis(2, las = 1, mgp = c(3, 0.75, 0), tcl = -0.45)
  mtext(y_axis_label, side = 2, line = 2.8)

  mtext(as.character(as.roman(seq_along(x_mids))),
        at = x_mids,
        side = 1,
        line = 0.3,
        cex = 0.9
        )
  mtext("Plate number", side = 1, line = 1.8)


  ## Indicate specific y axis ranges with colors
  color_scheme <- c(colorRampPalette(brewer.pal(9, "Greens"))(100)[[26]],
                    brewer.pal(9, "YlOrRd")[[2]],
                    colorRampPalette(brewer.pal(9, "Reds"))(100)[[21]]
                    )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = quality_ranges[[1]][[1]],
       ytop    = y_limits[[2]],
       col     = color_scheme[[1]],
       border  = NA
       )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = quality_ranges[[2]][[1]],
       ytop    = quality_ranges[[2]][[2]],
       col     = color_scheme[[2]],
       border  = NA
       )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = y_limits[[1]],
       ytop    = quality_ranges[[3]][[2]],
       col     = color_scheme[[3]],
       border  = NA
       )

  ## Draw horizontal and vertical indicator lines
  if (y_limits[[2]] != quality_ranges[[1]][[2]]) {
    abline(h = quality_ranges[[1]][[2]], col = "gray50", lty = "dotted")
  }
  abline(v = seq_len(length(x_mids) + 1) - 0.5,
         col = "gray70", lwd = 0.5
         )
  box()

  ## Plot the quality control metrics
  points(x   = x_positions,
         y   = data_vec,
         pch = 21,
         bg  = "black",
         cex = 0.7
         )

  par(old_mai)

  return(invisible(NULL))
}


GetQualityMetric <- function(input_df, UseFunction) {
  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  df_list <- split(input_df, plate_numbers_vec)
  rep1_vec <- vapply(df_list, UseFunction, use_column = "Raw_rep1", numeric(1))
  rep2_vec <- vapply(df_list, UseFunction, use_column = "Raw_rep2", numeric(1))
  results_mat <- cbind("rep1" = rep1_vec, "rep2" = rep2_vec)
  return(results_mat)
}


PlotZPrimes <- function(input_df) {
  z_primes_mat <- GetQualityMetric(input_df, Calculate_Z_Prime)
  PlotPlateQualities(z_primes_mat[, 1], z_primes_mat[, 2],
                     y_limits_include = c(0, 1, -0.2), y_axis_label = "Z' factor"
                     )
}


PlotSSMDControls <- function(input_df) {
  z_primes_mat <- GetQualityMetric(input_df, Calculate_SSMD_ctrls)
  PlotPlateQualities(z_primes_mat[, 1], z_primes_mat[, 2],
                     y_limits_include = c(0, 7),
                     y_axis_label = "SSMD (pos./neg. controls)",
                     quality_ranges = list(c(5, 7), c(3, 5), c(-Inf, 3))
                     )
}



ExportAllReplicateScatterPlots <- function(input_df) {

  use_dir <- file.path(output_dir, "Figures", "Replicate scatter plots")

  rep_columns <- grep("_rep", names(column_file_names), value = TRUE, fixed = TRUE)

  plot_height <- 4.5
  plot_ratio <- 0.45 / 0.21

  message("Exporting PDF plots...")

  pdf(file = file.path(use_dir, "Replicate scatter plots - flexible axes.pdf"),
      width = plot_height * plot_ratio, height = plot_height
      )
  for (use_column in rep_columns) {
    ReplicateScatter(input_df, rep1_column = use_column,
                     show_title = FormatPlotMath(long_column_labels[[use_column]]),
                     same_scale = FALSE
                     )
  }
  dev.off()


  pdf(file = file.path(use_dir, "Replicate scatter plots - fixed axes.pdf"),
      width = plot_height * plot_ratio, height = plot_height
      )
  for (use_column in rep_columns) {
    ReplicateScatter(input_df, rep1_column = use_column,
                     show_title = FormatPlotMath(long_column_labels[[use_column]]),
                     same_scale = TRUE
                     )
  }
  dev.off()

  message("Exporting PNG plots...")

  for (fixed_axes in c(FALSE, TRUE)) {
    for (i in seq_along(rep_columns)) {
      use_column <- rep_columns[[i]]
      file_name <- paste0("Replicate scatter plot - ", i,  ") ",
                          column_file_names[[use_column]], " - ",
                          if (fixed_axes) "fixed axes" else "flexible axes",
                          ".png"
                          )
      png(filename = file.path(use_dir, "Replicate scatter plots - PNGs", file_name),
          width = plot_height * plot_ratio, height = plot_height,
          units = "in", res = 600
          )
      ReplicateScatter(input_df, rep1_column = use_column,
                       show_title = FormatPlotMath(long_column_labels[[use_column]]),
                       same_scale = fixed_axes
                       )
      dev.off()
    }
  }
  return(invisible(NULL))
}


