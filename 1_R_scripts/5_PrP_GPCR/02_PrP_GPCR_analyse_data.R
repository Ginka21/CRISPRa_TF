# 2022-09-19


# Load packages and source code -------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
analysis_functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions", "2_Analyzing_data")
source(file.path(analysis_functions_dir, "01_Calculating_scores.R"))
source(file.path(analysis_functions_dir, "02_Processing_data.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP", "GPCRa")
output_dir <- file.path(project_dir, "4_output", "PrP", "Tables", "GPCRa")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "01_integrate_data.RData"))



# Avoid issues with taking the logarithm of negative values ---------------

are_low_rep1 <- PrP_df[, "Raw_rep1"] < 1
are_low_rep2 <- PrP_df[, "Raw_rep2"] < 1
are_gene <- !(is.na(PrP_df[, "Entrez_ID"]))
table(are_low_rep1)
table(are_low_rep2)
stopifnot(!(any(are_low_rep1 & (are_gene | PrP_df[, "Is_pos_ctrl"] | PrP_df[, "Is_NT_ctrl"]))))
stopifnot(!(any(are_low_rep2 & (are_gene | PrP_df[, "Is_pos_ctrl"] | PrP_df[, "Is_NT_ctrl"]))))
PrP_df[are_low_rep1, "Raw_rep1"] <- 1
PrP_df[are_low_rep2, "Raw_rep2"] <- 1




# Normalize by non-targeting controls -------------------------------------

PrP_df <- NormalizeWithNTControls(PrP_df, norm_method = "all NT")



# Calculate SSMD and derived statistics (p value, etc.) -------------------

PrP_df <- RunSSMDStats(PrP_df, norm_method = "all NT")




# Prepare hit list --------------------------------------------------------

hits_df_list <- CreateHitLists(PrP_df,
                               log2fc_column       = "Log2FC_rep1",
                               p_value_column      = "p_value_log2",
                               hit_strength_column = "Hit_strength_log2",
                               p_value_cutoff      = 0.05,
                               log2fc_cutoff       = log2(2)
                               )



# Examine criteria for defining hits --------------------------------------

use_df <- hits_df_list[["original_df"]]
are_gene <- !(is.na(use_df[, "Entrez_ID"]))

meet_p_val_cutoff  <- (use_df[, "p_value_log2"] < 0.05)
meet_log2fc_cutoff <- abs(use_df[, "Mean_log2FC"]) > log2(2)

meet_criteria <- meet_p_val_cutoff & meet_log2fc_cutoff
table(meet_criteria & are_gene)



# Check chosen cut-offs against the distribution of NT controls -----------

sum(meet_p_val_cutoff[use_df[, "Is_NT_ctrl"]])
sum(meet_log2fc_cutoff[use_df[, "Is_NT_ctrl"]])

range(use_df[, "p_value_log2"][use_df[, "Is_NT_ctrl"]])
NT_log2fc_range <- range(use_df[, "Mean_log2FC"][use_df[, "Is_NT_ctrl"]])
NT_log2fc_range
2^NT_log2fc_range

mean_NT <- mean(use_df[, "Mean_log2FC"][use_df[, "Is_NT_ctrl"]])
sd_NT <- sd(use_df[, "Mean_log2FC"][use_df[, "Is_NT_ctrl"]])
mean_NT + (c(-1, 1) * 3 * sd_NT)



# Export data -------------------------------------------------------------

write.csv(PrP_df,
          file = file.path(output_dir, "PrP_complete.csv"),
          row.names = FALSE, quote = FALSE
          )

exclude_columns <- c("Well_coords_384", grep("_96", names(PrP_df), fixed = TRUE, value = TRUE))
export_columns <- setdiff(names(hits_df_list[["reordered_df"]]),  exclude_columns)
write.csv(hits_df_list[["reordered_df"]][, export_columns],
          file = file.path(output_dir, "PrP_genes_and_NT_ordered.csv"),
          row.names = FALSE, quote = FALSE, na = ""
          )

write.csv(hits_df_list[["hits_df"]][, export_columns],
          file = file.path(output_dir,  "PrP_hits_only.csv"),
          row.names = FALSE, quote = FALSE, na = ""
          )




# Save data ---------------------------------------------------------------

save(PrP_df, file = file.path(r_data_dir, "02_analyse_data.RData"))


