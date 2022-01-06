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

metrics_in_order <- unique(columns_df[, "Metric"])
are_after <- seq_along(metrics_in_order) > which(metrics_in_order == "Fold-NT")
metrics_in_order <- unique(c(metrics_in_order[!(are_after)], "Log2FC", metrics_in_order[are_after]))

column_names_stripped <- sub("_Glo", "", columns_df[, "Column name"], fixed = FALSE)

new_order <- order(match(columns_df[, "Metric"], metrics_in_order),
                   match(column_names_stripped, column_names_stripped)
                   )
columns_df <- columns_df[new_order, ]
are_to_exclude <- grepl("_foldNT", columns_df[, "Column name"], fixed = TRUE) # p values calculated using foldNT values are not plausible
columns_df <- columns_df[!(are_to_exclude), ]
columns_df <- columns_df[columns_df[, "Metric"] != "t value", ]

column_labels <- columns_df[, "Label"]
names(column_labels) <- columns_df[, "Column name"]
names(column_labels) <- ifelse(columns_df[, "Replicates"] == "Yes",
                        paste0(names(column_labels), "_rep1"),
                        names(column_labels)
                        )


Glo_names <- c(
  "CellTiterGlo_raw"    = "Cell viability (CellTiterGlo values)",
  "CellTiterGlo_foldNT" = "Cell viability (normalized to NT controls)"
)

are_before <- seq_len(nrow(columns_df)) < which(columns_df[, "Column name"] == "SSMD_deltaNT")
column_labels <- c(column_labels[are_before], Glo_names, column_labels[!(are_before)])






