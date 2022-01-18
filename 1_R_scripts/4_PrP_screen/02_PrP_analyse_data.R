# 2022-01-18


# Load packages and source code -------------------------------------------

project_dir           <- "~/R_projects/CRISPRa_TF"
general_functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions", "2_Analyzing_data")
source(file.path(general_functions_dir, "01_Calculating_scores.R"))
source(file.path(general_functions_dir, "03_Processing_data.R"))



# Define folder path ------------------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir  <- file.path(project_dir,"4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "01_integrate_data.RData"))



# Avoid issues with taking the logarithm of negative values ---------------

are_low_rep1 <- PrP_df[, "Raw_rep1"] < 1
are_low_rep2 <- PrP_df[, "Raw_rep2"] < 1
are_gene <- !(is.na(PrP_df[, "Entrez_ID"]))
table(are_low_rep1)
table(are_low_rep2)
stopifnot(!(are_low_rep1 & (are_gene | PrP_df[, "Is_pos_ctrl"] | PrP_df[, "Is_NT_ctrl"])))
stopifnot(!(are_low_rep2 & (are_gene | PrP_df[, "Is_pos_ctrl"] | PrP_df[, "Is_NT_ctrl"])))
PrP_df[are_low_rep1, "Raw_rep1"] <- 1
PrP_df[are_low_rep2, "Raw_rep2"] <- 1



# Normalize by non-targeting controls -------------------------------------

PrP_df <- NormalizeWithNTControls(PrP_df)



# Calculate SSMD and derived statistics (p value, etc.) -------------------

PrP_df <- RunSSMDStats(PrP_df)



# Define hits -------------------------------------------------------------

are_gene <- !(is.na(PrP_df[, "Entrez_ID"]))
mean_logfc <- rowMeans(PrP_df[, c("Log2FC_rep1", "Log2FC_rep2")])

meet_p_val_cutoff  <- (PrP_df[, "p_value_log2"] < 0.05)
meet_log2fc_cutoff <- abs(mean_logfc) > log2(1.25)

meet_criteria <- meet_p_val_cutoff & meet_log2fc_cutoff



# Check chosen cut-offs against the distribution of NT controls -----------

# stopifnot(!(any(meet_criteria & PrP_df[, "Is_NT_ctrl"])))

sum(meet_p_val_cutoff[PrP_df[, "Is_NT_ctrl"]])
sum(meet_log2fc_cutoff[PrP_df[, "Is_NT_ctrl"]])

range(PrP_df[, "p_value_log2"][PrP_df[, "Is_NT_ctrl"]])
NT_log2fc_range <- range(mean_logfc[PrP_df[, "Is_NT_ctrl"]])
NT_log2fc_range
2^NT_log2fc_range

mean_NT <- mean(mean_logfc[PrP_df[, "Is_NT_ctrl"]])
sd_NT <- sd(mean_logfc[PrP_df[, "Is_NT_ctrl"]])
mean_NT + (c(-1, 1) * 3 * sd_NT)



# Prepare hit list --------------------------------------------------------

## Define additional columns that are useful for exporting the hit list
reordered_df <- PrP_df
reordered_df[, "Passes_cutoffs"]         <- meet_criteria
reordered_df[, "Mean_logFC"]             <- mean_logfc
reordered_df[, "p_value_log2_used"]      <- PrP_df[, "p_value_log2"]
reordered_df[, "Hit_strength_log2_used"] <- PrP_df[, "Hit_strength_log2"]


## Re-order columns to emphasize the data that was used for choosing hits,
## and re-order genes by their rank.
are_after <- seq_len(ncol(PrP_df)) > which(names(PrP_df) == "Is_pos_ctrl")
use_columns <- unique(c(names(PrP_df)[!(are_after)],
                        "Passes_cutoffs", "Mean_logFC", "p_value_log2_used",
                        "Hit_strength_log2_used",
                        names(PrP_df)[are_after]
                        ))

new_order <- order(meet_criteria,
                   abs(PrP_df[, "Hit_strength_log2"]),
                   decreasing = TRUE
                   )
reordered_df <- reordered_df[new_order, use_columns]

are_gene <- !(is.na(reordered_df[, "Entrez_ID"]))
are_selected <- are_gene | reordered_df[, "Is_NT_ctrl"]
reordered_df <- reordered_df[are_selected, ]
row.names(reordered_df) <- NULL


## Create a data frame containing only hit genes
are_hit <- (!(is.na(reordered_df[, "Entrez_ID"]))) & reordered_df[, "Passes_cutoffs"]
hits_df <- reordered_df[are_hit, ]
row.names(hits_df) <- NULL



# Export data -------------------------------------------------------------

write.csv(PrP_df,
          file = file.path(output_dir, "Tables", "PrP_complete.csv"),
          row.names = FALSE, quote = FALSE
          )

exclude_columns <- c("Well_coords_384", grep("_96", names(PrP_df), fixed = TRUE, value = TRUE))
write.csv(reordered_df[, !(names(reordered_df) %in% exclude_columns)],
          file = file.path(output_dir, "Tables", "PrP_genes_and_NT_ordered.csv"),
          row.names = FALSE, quote = FALSE, na = ""
          )

write.csv(hits_df,
          file = file.path(output_dir, "Tables", "PrP_hits_only.csv"),
          row.names = FALSE, quote = FALSE
          )


# Save data ---------------------------------------------------------------

save(PrP_df, file = file.path(r_data_dir, "02_analyse_data.RData"))


