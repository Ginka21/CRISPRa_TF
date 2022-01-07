# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_functions")
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


ScatterPlot <- function(x_vec,
                        y_vec,
                        use_limits = NULL,
                        point_size = 0.6,
                        label_y_axis = TRUE,
                        use_color = NULL,
                        top_label = NULL,
                        draw_regression = TRUE
                        ) {

  if (is.null(use_color)) {
    use_color <- "#000000"
  }

  stopifnot(length(x_vec) == length(y_vec))

  if (is.null(use_limits)) {
    use_limits = range(c(x_vec, y_vec))
  }

  use_limits <- use_limits + (diff(use_limits) * 0.04 * c(-1, 1))



  if (draw_regression) {
    ## Perform a linear regression
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


  use_mgp <- c(2.8, 0.55, 0)
  use_tcl <- -0.35
  plot(1,
       xlim = use_limits,
       ylim = use_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       type = "n"
       )
  axis(1, mgp = use_mgp, tcl = use_tcl, gap.axis = 0.15)
  mtext("Replicate 1", side = 1, line = 2.2)
  axis(2, las = 1, mgp = use_mgp, tcl = use_tcl)
  if (label_y_axis) {
    mtext("Replicate 2", side = 2, line = use_mgp[[1]])
  }
  if (!(is.null(top_label))) {
    mtext(top_label, line = 1.6, cex = par("cex"), font = 2)
  }

  mtext(VerticalAdjust(as.expression(corr_text)),
        line = 0.05, cex = par("cex"), font = 2
        )

  abline(a = 0, b = 1, col = "grey80", lty = "dotted")
  abline(h = 0, col = "gray90")
  abline(v = 0, col = "gray90")


  if (draw_regression) {
    polygon(c(new_df[, 1], rev(new_df[, 1])),
            c(conf_int_mat[, 2], rev(conf_int_mat[, 3])),
            col = Palify(use_color, fraction_pale = 0.8), border = NA
            )
    lines(new_df[, 1], conf_int_mat[, 1], col = use_color, lwd = 1.5)
  }

  box()

  points(x_vec,
         y_vec,
         pch = 16,
         col = adjustcolor(use_color, alpha.f = 0.5),
         cex = point_size * par("cex")
         )


  return(invisible(NULL))
}


MakeEmptyPlot <- function() {
  plot(1, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       type = "n", axes = FALSE, ann = FALSE
       )
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
  are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
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

  return(invisible(NULL))
}
#
#
# ReplicateScatter <- function(input_df,
#                              rep1_column,
#                              rep2_column = NULL,
#                              show_title  = "",
#                              point_size = 0.6
#                              ) {
#
#   if (is.null(rep2_column)) {
#     rep2_column <- sub("rep1", "rep2", rep1_column, fixed = TRUE)
#   }
#
#   are_genes <- !(is.na(input_df[, "Entrez_ID"]))
#
#   use_limits <- range(input_df[are_genes, c(rep1_column, rep2_column)])
#
#
#   cor_results <- cor.test(input_df[are_genes, rep1_column],
#                           input_df[are_genes, rep2_column]
#                           )
#   r_value <- cor_results[["estimate"]][[1]]
#
#   old_mai <- par(mai = c(0.9, 1, 0.6, 1.4))
#
#   plot(GBA_df[are_genes, rep1_column],
#        GBA_df[are_genes, rep2_column],
#        xlim = use_limits,
#        ylim = use_limits,
#        xlab = "Replicate 1",
#        ylab = "Replicate 2",
#        main = show_title,
#        cex.main = 1.1,
#        font.main = 1,
#        pch  = 16,
#        cex  = point_size,
#        col  = adjustcolor("black", alpha.f = 0.5),
#        las  = 1,
#        mgp  = c(2.9, 0.65, 0),
#        tcl  = -0.45
#        )
#
#   x_start <- 1.03
#
#   text(x      = grconvertX(x_start, from = "npc", to = "user"),
#        y      = grconvertY(0.9, from = "npc", to = "user"),
#        labels = expression(bold("Pearson's " * bolditalic("r"))),
#        adj    = c(0, 0),
#        xpd    = NA
#        )
#
#   text(x      = grconvertX(x_start, from = "npc", to = "user"),
#        y      = grconvertY(0.82, from = "npc", to = "user"),
#        labels = bquote(bold("= " * .(as.character(round(r_value, digits = 2))))),
#        adj    = c(0, 0),
#        xpd    = NA
#        )
#
#   par(old_mai)
#
#   return(invisible(NULL))
# }
#



# Calculate correlation between replicates --------------------------------

ReplicateScatter(GBA_df, "Raw_rep1")
ReplicateScatter(GBA_df, "PercActivation_log2_Glo_rep1")


rep_columns <- grep("_rep", names(column_labels), value = TRUE, fixed = TRUE)

plot_height <- 4.5
plot_ratio <- 0.45 / 0.21


pdf(file = file.path(output_dir, "Figures", "Replicate scatter plots", "Replicate scatter plots - flexible axes.pdf"),
    width = plot_height * plot_ratio, height = plot_height
    )
for (use_column in rep_columns) {
  ReplicateScatter(GBA_df, rep1_column = use_column,
                   show_title = column_labels[[use_column]], same_scale = FALSE
                   )
}
dev.off()


pdf(file = file.path(output_dir, "Figures", "Replicate scatter plots", "Replicate scatter plots - fixed axes.pdf"),
    width = plot_height * plot_ratio, height = plot_height
    )
for (use_column in rep_columns) {
  ReplicateScatter(GBA_df, rep1_column = use_column,
                   show_title = column_labels[[use_column]], same_scale = TRUE
                   )
}
dev.off()


for (fixed_axes in c(FALSE, TRUE)) {
  for (i in seq_along(rep_columns)) {
    use_column <- rep_columns[[i]]
    file_name <- paste0("Replicate scatter plot - ", i,  ") ",
                        sub("_rep1", "", use_column, fixed = TRUE), " - ",
                        if (fixed_axes) "fixed axes" else "flexible axes",
                        ".png"
                        )
    png(file = file.path(output_dir, "Figures", "Replicate scatter plots", "Replicate scatter plots - PNGs", file_name),
        width = plot_height * plot_ratio, height = plot_height,
        units = "in", res = 600
        )
    ReplicateScatter(GBA_df, rep1_column = use_column,
                     show_title = column_labels[[use_column]],
                     same_scale = fixed_axes
                     )
    dev.off()
  }
}



# Calculate Z' Factor for every plate ------------------------------------

z_prime_vec_1 <- vapply(split_df_list,
                       Calculate_Z_Prime,
                       use_column = "Raw_rep1",
                       numeric(1)
                       )

z_prime_vec_2 <- vapply(split_df_list,
                       Calculate_Z_Prime,
                       use_column = "Raw_rep2",
                       numeric(1)
                       )



# Calculate SSMD for every plate (using controls) -------------------------

SSMD_vec_1 <- vapply(split_df_list,
                       Calculate_SSMD_ctrls,
                       use_column = "Raw_rep1",
                       numeric(1)
                       )

SSMD_vec_2 <- vapply(split_df_list,
                       Calculate_SSMD_ctrls,
                       use_column = "Raw_rep2",
                       numeric(1)
                       )


