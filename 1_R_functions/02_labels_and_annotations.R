# 2021-12-29


# Load packages and source code -------------------------------------------

library("readxl")



# Define folder path ------------------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
input_dir   <- file.path(project_dir, "2_input")




# Read in data ------------------------------------------------------------

columns_df <- data.frame(read_excel(file.path(input_dir, "All_metrics.xlsx")),
                         check.names = FALSE, stringsAsFactors = FALSE
                         )




# Define labels -----------------------------------------------------------

column_labels <- columns_df[, "Label"]

metrics_in_order <- unique(columns_df[, "Metric"])
are_after <- seq_along(metrics_in_order) > which(metrics_in_order == "Fold-NT")
metrics_in_order <- unique(c(metrics_in_order[!(are_after)], "Log2FC", metrics_in_order[are_after]))

new_order <- order(match(columns_df[, "Metric"], metrics_in_order))
columns_df <- columns_df[new_order, ]

column_labels <- ifelse(columns_df[, "Replicates"] == "Yes",
                        paste0(column_labels, "_rep1"),
                        column_labels
                        )
names(column_labels) <- columns_df[, "Column name"]


Glo_names <- c(
  "CellTiterGlo_raw"    = "Cell viability (CellTiterGlo values)",
  "CellTiterGlo_foldNT" = "Cell viability (normalized to NT controls)"
)

are_Glo <- columns_df[, "Standardized by viability (CellTiterGlo)"] %in% "Yes"

column_labels <- c(column_labels[!(are_Glo)], Glo_names, column_labels[are_Glo])


