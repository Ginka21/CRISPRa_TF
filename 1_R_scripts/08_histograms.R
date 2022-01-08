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


# Define Functions --------------------------------------------------------

DrawHist <- function (input_df,
                      use_column,
                      breaks = use_breaks,
                      x_axis_label = NULL,
                      color = NA
                      ){

  if (is.null(x_axis_label)) {
      if (("column_labels" %in% ls(envir = globalenv())) &&
          (use_column %in% names(column_labels))
          ) {
         x_axis_label <- column_labels[[use_column]]
      } else {
         x_axis_label <- use_column
      }
   }

   are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
   are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
   are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
   are_valid   <- are_NT | are_posctrl | are_gene

   has_replicates <- grepl("_rep", use_column, fixed = TRUE)
   if (has_replicates) {
     rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
     combined_vec <- c(input_df[, use_column], input_df[, rep2_column])

     gene_vec     <- combined_vec[rep(are_gene, 2)]
     pos_ctrl_vec <- combined_vec[rep(are_posctrl, 2)]
     NT_ctrl_vec  <- combined_vec[rep(are_NT, 2)]
     valid_vec    <- combined_vec[rep(are_valid, 2)]

   } else {
     numeric_vec  <- input_df[, use_column]
     gene_vec     <- numeric_vec[are_gene]
     pos_ctrl_vec <- numeric_vec[are_pos_ctrl]
     NT_ctrl_vec  <- numeric_vec[are_NT_ctrl]
     valid_vec    <- numeric_vec[are_valid]
   }

   pos_ctrl_color   <- brewer.pal(5, "Reds")[[4]]
   NT_ctrl_color    <- brewer.pal(5, "Blues")[[3]]
   gene_color       <- brewer.pal(5, "Greens")[[3]]

   use_margin <- c(7, 4.6, 3.8, 7.5)
   old_mar    <- par(mar = use_margin)
   use_mgp    <- c(3, 0.65, 0)

   # If the minimum and maximum are closer together than the distance
   # from the minimum to zero, then the x axis should not include zero
   far_from_zero <- diff(range(valid_vec)) < min(valid_vec)
   if (far_from_zero) {
      x_limits <- range(valid_vec)
   } else {
      x_limits <- range(c(0, valid_vec))
   }
   x_space <- diff(x_limits) * 0.04
   x_limits <- x_limits + (x_space * c(-1, 1))

   # If the minimum value is not negative or close to zero, the x range should stop at zero
   if (!(far_from_zero) && (min(valid_vec) >= x_space)) {
      x_limits[[1]] <- 0
   }

   hist_results <- hist(1,
                       breaks = 1000,
                       color  = use_color,
                       xaxs   = "i",
                       yaxs   = "i",
                       border = NA,
                       main   = "",
                       xlab   = x_axis_label,
                       mgp    = c(2.5, 0.5, 0),
                       freq   = TRUE,
                       las    = 1,
                       ylab   = "Frequency",
                       type   = "n"
                       )

                       hist(gene_vec,
                       col    = gene_vec_color,
                       breaks = use_breaks,
                       add    = TRUE
                       )

                       hist(pos_ctrl_vec,
                       col    = pos_ctrl_color,
                       breaks = use_breaks,
                       add    = TRUE
                       )

                       hist(neg_ctrl_vec,
                       col    = neg_ctrl_color,
                       breaks = use_breaks,
                       add    = TRUE
                       )

# Add legend
  controls_labels <- list(
      "Gene" = c("Genes in ", "CRISPRa", "library"),
      "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
      "NT"   = c("Non-targeting", "controls")
   )


   DrawSideLegend(labels_list = controls_labels,
                  use_pch = 15,
                  use_colors = c(gene_color, pos_ctrl_color, NT_ctrl_color)
                  )

  return(invisible(hist_results))
}


# my_hist <- function(numeric_vec, ..., width = 1000, anchor = min(numeric_vec), range_lables = FALSE) {
#     offset <- anchor %% width
#     low <- width * (min(numeric_vec) - offset) %/% width + offset
#     high <- width * (max(numeric_vec) + width - offset - 1) %/% width + offset
#     br <-  seq(low, high, width)
#     hist(numeric_vec, breaks=br, ..., xaxt = "n")
#     if (range_lables) {
#       labs <- paste(br[seq_along(br[-1])] + c(0, rep(1, length(br) - 2)), br[-1], sep = "-")
#       at <- br[seq_along(br[-1])] + width / 2
#       axis(1, at, labs, las = 2)
#     } else {
#       axis(side=1, at = br)
#   }
# }

# Create Histogram for Viability -----------------------------------------

input_df     <- GBA_df
are_NT       <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
are_pos_ctrl <- input_df[, "Target_flag"] %in% "Pos. control"
are_gene     <- !(is.na(input_df[, "Entrez_ID"]))
are_noempty  <- are_NT | are_pos_ctrl | are_gene


lum_gene_vec     <- input_df[are_gene, "Luminescence"]
lum_pos_ctrl_vec <- input_df[are_pos_ctrl, "Luminescence"]
lum_NT_ctrl_vec  <- input_df[are_NT, "Luminescence"]
lum_noempty_vec  <- input_df[are_noempty, "Luminescence"]


pos_ctrl_color   <- brewer.pal(5, "Reds")[[4]]
NT_ctrl_color    <- brewer.pal(5, "Blues")[[3]]
gene_color       <- brewer.pal(5, "Greens")[[3]]

use_margin <- c(7, 4.6, 3.8, 7.5)
old_mar    <- par(mar = use_margin)
use_mgp    <- c(3, 0.65, 0)
x_limits   <- range(lum_noempty_vec)
x_limits   <- x_limits + (diff(x_limits) * 0.04 * c(-1, 1))
use_breaks <- seq(min(lum_noempty_vec), max(lum_noempty_vec), by=((max(lum_noempty_vec) - min(lum_noempty_vec))/(length(lum_noempty_vec))))

# Plot First distribution
hist(lum_gene_vec,
     breaks = use_breaks,
     xlim   = x_limits,
     col    = adjustcolor(gene_color, alpha.f = 0.8),
     border = NA,
     las    = 1,
     xlab   = "Absolute luminescent values",
     ylab   = "Frequency",
     main   = "Cell Viability (CellTiterGlo)" )

# Second with add=T to plot on top
hist(lum_pos_ctrl_vec,
     breaks = use_breaks,
     xlim   = x_limits,
     col    = adjustcolor(pos_ctrl_color, alpha.f = 0.8),
     border = NA,
     add    = TRUE
     )

# Third with add=T to plot on top
hist(lum_NT_ctrl_vec,
     breaks = use_breaks,
     xlim   = x_limits,
     col    = adjustcolor(NT_ctrl_color, alpha.f = 0.8),
     border = NA,
     add    = TRUE
     )

# Add legend
  controls_labels <- list(
      "Gene" = c("Genes in ", "CRISPRa", "library"),
      "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
      "NT"   = c("Non-targeting", "controls")
   )



   DrawSideLegend(labels_list = controls_labels,
                  use_pch = 15,
                  use_colors = c(gene_color, pos_ctrl_color, NT_ctrl_color)
                  )


# Create Histograms for controls and samples ------------------------------

are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
are_noempty <- are_NT | are_posctrl | are_gene

combined_vec <- c(input_df[, "GBA_rep1_absolute"], input_df[, "GBA_rep2_absolute"])

gene_vec     <- combined_vec[rep(are_gene, 2)]
pos_ctrl_vec <- combined_vec[rep(are_posctrl, 2)]
NT_ctrl_vec  <- combined_vec[rep(are_NT, 2)]
noempty_vec  <- combined_vec[rep(are_noempty, 2)]

use_margin <- c(7, 4.6, 3.8, 7.5)
old_mar    <- par(mar = use_margin)
use_mgp    <- c(3, 0.65, 0)
x_range    <- c(0, max(noempty_vec) * 1.05)

pos_ctrl_color   <- brewer.pal(9, "Reds")[[7]]
NT_ctrl_color    <- brewer.pal(9, "Blues")[[7]]
gene_color       <- brewer.pal(9, "Greens")[[7]]

# usebreaks <- as.matrix(noempty_vec) [["breaks"]]
usebreaks <- 50

# Plot First distribution
hist(gene_vec,
     breaks = usebreaks,
     xlim   = x_range,
     col    = adjustcolor(gene_color, alpha.f = 0.8),
     border = NA,
     las    = 1,
     xlab   = "Absolute Fluorescent values",
     ylab   = "Frequency",
     main   = "Absolute Fluorecence"
     )

# add = TRUE to plot on top
hist(pos_ctrl_vec,
     breaks = usebreaks,
     xlim   = c(0, max_xlim),
     col    = adjustcolor(pos_ctrl_color, alpha.f = 0.8),
     border = NA,
     add    = TRUE)


hist(NT_ctrl_vec,
     breaks = usebreaks,
     xlim   = c(0, max_xlim),
     col    = adjustcolor(NT_ctrl_color, alpha.f = 0.8),
     border = NA,
     add    = TRUE
     )

# Add legend
  controls_labels <- list(
      "Gene" = c("Genes in ", "CRISPRa", "library"),
      "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
      "NT"   = c("Non-targeting", "controls")
   )


   DrawSideLegend(labels_list = controls_labels,
                  use_pch = 15,
                  use_colors = c(gene_color, pos_ctrl_color, NT_ctrl_color)
                  )


# Draw histograms ---------------------------------------------------------

# categ_mat <- sapply(levels(essential_df[["Four_categories"]]), function(x) {
#   are_this_category <- essential_df[["Four_categories"]] %in% x
#   CRISPR_effects_df[["Entrez_ID"]] %in% essential_df[["Entrez_ID"]][are_this_category]
# })
#
#
#
# DrawHistogram <- function(numeric_vec,
#                           add = FALSE,
#                           hist_color = brewer.pal(9, "Blues")[[8]],
#                           use_breaks = 1000
#                           ) {
#   points_alpha <- 0.5
#   alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
#   if (!(is.na(hist_color))) {
#     use_color <- paste0(hist_color, alpha_hex)
#   } else {
#     use_color <- NA
#   }
#   hist_results <- hist(numeric_vec,
#                        breaks = use_breaks,
#                        col    = use_color,
#                        border = NA,
#                        main   = "Depmap \u2013 all cell lines",
#                        xlab   = "CRISPR knockout fitness effect",
#                        mgp    = c(2.5, 0.5, 0),
#                        freq   = TRUE,
#                        add    = add,
#                        axes   = FALSE,
#                        ylab   = ""
#                        )
#   if (!(add)) {
#     box(bty = "l")
#     x_axis_pos <- pretty(par("usr")[c(1, 2)], n = 10)
#     axis(1, at = x_axis_pos, mgp = c(2.5, 0.55, 0), tcl = -0.45)
#   }
#   return(invisible(hist_results))
# }
#
#
#
# for (make_PNG in c(TRUE, FALSE)) {
#
#   if (make_PNG) {
#     png(file = file.path(file_output_directory, "Histograms - gene effects.png"),
#         height = 6, width = 8, units = "in", res = 900
#         )
#   }
#
#   hist_breaks <- DrawHistogram(as.matrix(CRISPR_effects_df[, 4:ncol(CRISPR_effects_df)]),
#                                hist_color = NA
#                                )[["breaks"]]
#
#   abline(v = seq(-0.1, 0.1, by = 0.1), col = c("gray75", "gray50"), lty = "dashed")
#
#
#   DrawHistogram(as.matrix(CRISPR_effects_df[, 4:ncol(CRISPR_effects_df)]),
#                 hist_color = brewer.pal(9, "Greys")[[4]],
#                 add = TRUE, use_breaks = hist_breaks
#                 )
#
#
#   DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Non-essential"], 4:ncol(CRISPR_effects_df)]),
#                 add = TRUE, hist_color = brewer.pal(9, "Greens")[[8]], use_breaks = hist_breaks
#                 )
#
#   DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Intermediate"], 4:ncol(CRISPR_effects_df)]),
#                 add = TRUE, hist_color = "#aa6c39", use_breaks = hist_breaks
#                 )
#
#   DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Essential"], 4:ncol(CRISPR_effects_df)]),
#                 add = TRUE, hist_color = brewer.pal(9, "Reds")[[8]], use_breaks = hist_breaks
#                 )
#
#
#   legend("topleft",
#          legend     = c("All genes", "Non-essential", "Intermediate", "Essential"),
#          fill       = c(brewer.pal(9, "Greys")[[4]],
#                         brewer.pal(9, "Greens")[[8]],
#                         "#aa6c39",
#                         brewer.pal(9, "Reds")[[8]]
#                         ),
#          border    = NA,
#          bty       = "n",
#          y.intersp = 1.1
#          )
#
#   if (make_PNG) {
#     dev.off()
#   }
#
# }
#
