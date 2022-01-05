# 2021-12-27


# Load packages and source code -------------------------------------------

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_functions")
source(file.path(functions_dir, "01_calculating_scores.R"))



# Define folder path ------------------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_integrate_data.RData"))



# Normalize by non-targeting controls -------------------------------------
##Follow the split-apply-combine strategy

plate_numbers_vec <- as.integer(as.roman(GBA_df[, "Plate_number_384"]))
split_df_list <- split(GBA_df, plate_numbers_vec)

split_df_list <- lapply(split_df_list, function(x) {

  are_NT <- x[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  are_pos <- x[, "Target_flag"] %in% "Pos. control"

  # Normalize fluorescent values
  dup1_median_NT  <- median(x[are_NT, "Raw_rep1"])
  dup2_median_NT  <- median(x[are_NT, "Raw_rep2"])
  dup1_median_pos <- median(x[are_pos, "Raw_rep1"])
  dup2_median_pos <- median(x[are_pos, "Raw_rep2"])

  x[, "FoldNT_rep1"]         <- x[, "Raw_rep1"] / dup1_median_NT
  x[, "FoldNT_rep2"]         <- x[, "Raw_rep2"] / dup2_median_NT
  x[, "DeltaNT_rep1"]        <- x[, "Raw_rep1"] - dup1_median_NT
  x[, "DeltaNT_rep2"]        <- x[, "Raw_rep2"] - dup2_median_NT
  x[, "PercActivation_rep1"] <- x[, "DeltaNT_rep1"] / (dup1_median_pos - dup1_median_NT)
  x[, "PercActivation_rep2"] <- x[, "DeltaNT_rep2"] / (dup2_median_pos - dup2_median_NT)


  # # Log transform fluorescent values
  x[, "Raw_log2_rep1"] <- log2(x[, "Raw_rep1"])
  x[, "Raw_log2_rep2"] <- log2(x[, "Raw_rep2"])
  dup1_NT_vec_log2 <- x[are_NT, "Raw_log2_rep1"]
  dup2_NT_vec_log2 <- x[are_NT, "Raw_log2_rep2"]
  x[, "Log2FC_rep1"] <- x[, "Raw_log2_rep1"] - median(dup1_NT_vec_log2)
  x[, "Log2FC_rep2"] <- x[, "Raw_log2_rep2"] - median(dup2_NT_vec_log2)




  # Normalize luminescent values
  lum_NT_vec                 <- x[are_NT, "CellTiterGlo_raw"]
  x[, "CellTiterGlo_foldNT"] <- x[, "CellTiterGlo_raw"] / median(lum_NT_vec)

  # Standardized Fluorescence by Luminescence
  x[, "Raw_Glo_rep1"] <- x[, "Raw_rep1"] /
                                      x[, "CellTiterGlo_foldNT"]

  x[, "Raw_Glo_rep2"] <- x[, "Raw_rep2"] /
                                      x[, "CellTiterGlo_foldNT"]


  # Log transform fluorescent values standardized by luminescence
  x[, "Raw_log2_Glo_rep1"] <- log2(x[, ("Raw_Glo_rep1")])
  x[, "Raw_log2_Glo_rep2"] <- log2(x[, ("Raw_Glo_rep2")])

  # Normalize Glo standardized values
  dup1_median_NT <- median(x[are_NT, "Raw_log2_Glo_rep1"])
  dup2_median_NT <- median(x[are_NT, "Raw_log2_Glo_rep2"])
  x[, "GBA_rep1_norm_Glo"] <- x[, "GBA_rep1_Glo_standardized"] / dup1_median_NT
  x[, "GBA_rep2_norm_Glo"] <- x[, "GBA_rep2_Glo_standardized"] / dup2_median_NT
  x[, "GBA_rep1_diff_Glo"] <- x[, "GBA_rep1_Glo_standardized"] - dup1_median_NT
  x[, "GBA_rep2_diff_Glo"] <- x[, "GBA_rep2_Glo_standardized"] - dup2_median_NT

  return(x)
})

GBA_df <- do.call(rbind.data.frame, c(split_df_list, make.row.names = FALSE))



stopifnot(identical(Calculate_SSMD_var(GBA_df, "Raw_rep1"),
                    Calculate_SSMD_var_old(GBA_df, "Raw_rep1")
                    ))

stopifnot(identical(Calculate_SSMD_var(GBA_df, "Raw_rep1", log2FC = TRUE),
                    Calculate_SSMD_var_old(GBA_df, "Raw_rep1", log2FC = TRUE)
                    ))




# Calculate log-fold change -----------------------------------------------

GBA_df[, "GBA_rep1_logFC"]     <- log2(GBA_df[, "GBA_rep1_normalized"])
GBA_df[, "GBA_rep2_logFC"]     <- log2(GBA_df[, "GBA_rep2_normalized"])
GBA_df[, "GBA_rep1_Glo_logFC"] <- log2(GBA_df[, "GBA_rep1_norm_Glo"])
GBA_df[, "GBA_rep2_Glo_logFC"] <- log2(GBA_df[, "GBA_rep2_norm_Glo"])


# Calculate SSMD ----------------------------------------------------------

GBA_df[, "SSMD_MM_paired"]     <- Calculate_SSMD_var(GBA_df, "Raw_rep1")
GBA_df[, "SSMD_MM_paired_Glo"] <- Calculate_SSMD_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "SSMD_log2"]          <- Calculate_SSMD_var(GBA_df, "Raw_rep1", log2FC = TRUE)
GBA_df[, "SSMD_log2_Glo"]      <- Calculate_SSMD_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)

GBA_df[, "SSMD_test"]           <- Calculate_SSMD_var2(GBA_df, "Raw_rep1")
GBA_df[, "SSMD_test_Glo"]       <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "SSMD_log2_test"]      <- Calculate_SSMD_var2(GBA_df, "Raw_rep1", log2FC = TRUE)
GBA_df[, "SSMD_log2_Glo_test"]  <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)
GBA_df[, "SSMD_test_log2_diff"] <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_log2_diff")

# Calculate T Score ----------------------------------------------------------

GBA_df[, "T_value"]         <- Calculate_T_var(GBA_df, "Raw_rep1")
GBA_df[, "T_value_Glo"]     <- Calculate_T_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "T_value_log"]     <- Calculate_T_var(GBA_df, "Raw_rep1", log2FC = TRUE)
GBA_df[, "T_value_Glo_log"] <- Calculate_T_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)


# Calculate P value ----------------------------------------------------------

GBA_df[, "P_value"]         <- Calculate_P_var(GBA_df, "Raw_rep1")
GBA_df[, "P_value_Glo"]     <- Calculate_P_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "P_value_log"]     <- Calculate_P_var(GBA_df, "Raw_rep1", log2FC = TRUE)
GBA_df[, "P_value_Glo_log"] <- Calculate_P_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)



# Calculate hit strength --------------------------------------------------

GBA_df[, "Hit_strength"]     <- rowMeans(log2(GBA_df[, c("GBA_rep1_normalized", "GBA_rep2_normalized")])) *
                                -log10(GBA_df[, "P_value_log"])
GBA_df[, "Hit_strength_Glo"] <- rowMeans(log2((GBA_df[, c("GBA_rep1_norm_Glo", "GBA_rep2_norm_Glo")]))) *
                                -log10(GBA_df[, "P_value_Glo_log"])




# Export Data -------------------------------------------------------------

write.csv(GBA_df,
          file = file.path(output_dir, "Tables", "GBA_complete.csv"),
          row.names = FALSE, quote = FALSE
          )



# Save data ---------------------------------------------------------------

save(GBA_df, file = file.path(r_data_dir, "03_analyse_data.RData"))
