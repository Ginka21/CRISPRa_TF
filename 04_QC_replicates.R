# 2021-12-27


# Load packages and source code -------------------------------------------

project_dir   <- "~/R_projects/Screen"
functions_dir <- file.path(project_dir, "1_R_functions")
source(file.path(functions_dir, "01_calculating_scores.R"))
source(file.path(functions_dir, "02_labels_and_annotations.R"))
library("RColorBrewer")


# Define folder path ------------------------------------------------------

input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")


# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "03_analyse_data.RData"))
plate_numbers_vec <- as.integer(as.roman(GBA_df[, "Plate_number_384"]))
split_df_list <- split(GBA_df, plate_numbers_vec)


# Define Functions --------------------------------------------------------

ReplicateScatter <- function(input_df,
                             rep1_column,
                             rep2_column = NULL,
                             show_title  = ""
                             ) {

  if (is.null(rep2_column)) {
    rep2_column <- sub("rep1", "rep2", rep1_column, fixed = TRUE)
  }

  are_genes <- !(is.na(input_df[, "Entrez_ID"]))

  use_limits <- range(input_df[are_genes, c(rep1_column, rep2_column)])


  cor_results <- cor.test(input_df[are_genes, rep1_column],
                          input_df[are_genes, rep2_column]
                          )
  r_value <- cor_results[["estimate"]][[1]]

  old_mai <- par(mai = c(0.9, 1, 0.6, 1.4))

  plot(GBA_df[are_genes, rep1_column],
       GBA_df[are_genes, rep2_column],
       xlim = use_limits,
       ylim = use_limits,
       xlab = "Replicate 1",
       ylab = "Replicate 2",
       main = show_title,
       cex.main = 1.1,
       font.main = 1,
       pch  = 16,
       col  = adjustcolor("black", alpha.f = 0.5),
       las  = 1,
       mgp  = c(2.9, 0.65, 0),
       tcl  = -0.45
       )

  x_start <- 1.03

  text(x      = grconvertX(x_start, from = "npc", to = "user"),
       y      = grconvertY(0.9, from = "npc", to = "user"),
       labels = expression(bold("Pearson's " * bolditalic("r"))),
       adj    = c(0, 0),
       xpd    = NA
       )

  text(x      = grconvertX(x_start, from = "npc", to = "user"),
       y      = grconvertY(0.82, from = "npc", to = "user"),
       labels = bquote(bold("= " * .(as.character(round(r_value, digits = 2))))),
       adj    = c(0, 0),
       xpd    = NA
       )

  par(old_mai)

  return(invisible(NULL))
}




# Calculate correlation between replicates --------------------------------

ReplicateScatter(GBA_df, "GBA_rep1_absolute")


rep_columns <- grep("_rep", names(column_labels), value = TRUE, fixed = TRUE)


pdf(file = file.path(output_dir, "Figures", "Replicate scatter plots.pdf"),
    width = 6, height = 5.1
    )
for (use_column in rep_columns) {
  ReplicateScatter(GBA_df, rep1_column = use_column, show_title = column_labels[[use_column]])
}
dev.off()



for (i in seq_along(rep_columns)) {
  use_column <- rep_columns[[i]]
  file_name <- paste0("Replicate scatter plot - ", i,  ") ",
                      sub("_rep1", "", use_column, fixed = TRUE),
                      ".png"
                      )
  png(file = file.path(output_dir, "Figures", "Replicate scatter plots - PNGs", file_name),
      width = 6, height = 5.1, units = "in", res = 600
      )
  ReplicateScatter(GBA_df, rep1_column = use_column, show_title = column_labels[[use_column]])
  dev.off()
}


# Calculate Z' Factor  for every plate ------------------------------------

z_prime_vec_1 <- vapply(split_df_list,
                       Calculate_Z_Prime,
                       use_column = "GBA_rep1_absolute",
                       numeric(1)
                       )

z_prime_vec_2 <- vapply(split_df_list,
                       Calculate_Z_Prime,
                       use_column = "GBA_rep2_absolute",
                       numeric(1)
                       )


stop()

# Create Histogram for Viability -----------------------------------------

input_df     <- GBA_df
are_NT       <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
are_pos_ctrl <- input_df[, "Target_flag"] %in% "Pos. control"
are_gene     <- !(is.na(input_df[, "Entrez_ID"]))
are_noempty  <- are_NT | are_pos_ctrl | are_gene

pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
NT_ctrl_color  <- brewer.pal(5, "Blues")[[3]]


lum_gene_vec     <- input_df[are_gene, "Luminescence"]
lum_pos_ctrl_vec <- input_df[are_pos_ctrl, "Luminescence"]
lum_NT_ctrl_vec  <- input_df[are_NT, "Luminescence"]
lum_noempty_vec  <- input_df[are_noempty, "Luminescence"]

max_xlim         <- max(lum_noempty_vec) + 2000
pos_ctrl_color   <- brewer.pal(5, "Reds")[[4]]
NT_ctrl_color    <- brewer.pal(5, "Blues")[[3]]
gene_color       <- brewer.pal(5, "Greens")[[3]]

# Plot First distribution
hist(lum_gene_vec,
     breaks = 30,
     xlim   = c(0, max_xlim),
     col    = gene_color,
     xlab   = "Absolute luminescent values",
     ylab   = "Frequency",
     main   = "Cell Viability (CellTiterGlo)" )

# Second with add=T to plot on top
hist(lum_pos_ctrl_vec,
     breaks = 30,
     xlim   = c(0, max_xlim),
     col    = pos_ctrl_color,
     add    = T
     )

# Third with add=T to plot on top
hist(lum_NT_ctrl_vec,
     breaks = 30,
     xlim = c(0, max_xlim),
     col = NT_ctrl_color,
     add = T
     )

# Add legend
legend("topleft",
       legend=c("Genes CRISPRa library","Pos. ctrl", "NT ctrl"),
       col=c(gene_color, pos_ctrl_color, col = NT_ctrl_color),
       pt.cex=2,
       pch=15 )

# Create Histograms for controls and samples ------------------------------

are_NT      <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
are_posctrl <- input_df[, "Target_flag"] %in% "Pos. control"
are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
are_noempty <- are_NT | are_posctrl | are_gene

rep1_df <- input_df[, "GBA_rep1_absolute"]
rep2_df <- input_df[, "GBA_rep2_absolute"]
combined_vec <- lapply(rep1_df, rep2_df, SIMPLIFY = FALSE)

gene_vec     <- combined_vec[are_gene]
pos_ctrl_vec <- combined_vec[are_posctrl]
NT_ctrl_vec  <- combined_vec[are_NT]
noempty_vec  <- combined_vec[are_noempty]


max_xlim         <- range(0, max(noempty_vec) + 500)
max_ylim         <- length(noempty_vec)
pos_ctrl_color   <- brewer.pal(5, "Reds")[[4]]
NT_ctrl_color    <- brewer.pal(5, "Blues")[[3]]
gene_color       <- brewer.pal(5, "Greens")[[3]]
custombreaks     <- 17

# Plot First distribution
hist(gene_vec,
     breaks = custombreaks,
     xlim = c(0, max_xlim),
     # ylim = max_ylim,
     col = gene_color,
     xlab = "Absolute Fluorescent values",
     ylab ="Frequency",
     main ="Absolute Fluorecence" )

# Second with add=T to plot on top
hist(pos_ctrl_vec,
     breaks = custombreaks,
     xlim= c(0, max_xlim),
     col = pos_ctrl_color,
     add=T)

# Third with add=T to plot on top
hist(NT_ctrl_vec,
     breaks = 5,
     xlim=c(0, max_xlim),
     col = NT_ctrl_color,
     add=T)

# Add legend
legend("topright", legend=c("Genes CRISPRa library","Pos. ctrl", "NT ctrl"), col=c(gene_color,
      pos_ctrl_color, col = NT_ctrl_color), pt.cex=2, pch=15 )

