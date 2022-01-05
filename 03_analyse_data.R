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

  # Normalize fluorescent values
  dup1_NT_vec       <- x[are_NT, "GBA_rep1_absolute"]
  dup2_NT_vec       <- x[are_NT, "GBA_rep2_absolute"]
  dup1_median_NT    <- median(dup1_NT_vec)
  dup2_median_NT    <- median(dup2_NT_vec)

  x[, "GBA_rep1_normalized"] <- x[, "GBA_rep1_absolute"] / dup1_median_NT
  x[, "GBA_rep2_normalized"] <- x[, "GBA_rep2_absolute"] / dup2_median_NT
  x[, "GBA_rep1_diff"]       <- x[, "GBA_rep1_absolute"] - dup1_median_NT
  x[, "GBA_rep2_diff"]       <- x[, "GBA_rep2_absolute"] - dup2_median_NT

  # # Log transform fluorescent values
  x[, "GBA_rep1_log2"]   <- log2(x[, "GBA_rep1_absolute"])
  x[, "GBA_rep2_log2"]   <- log2(x[, "GBA_rep2_absolute"])
  dup1_NT_vec_log2       <- x[are_NT, "GBA_rep1_log2"]
  dup2_NT_vec_log2       <- x[are_NT, "GBA_rep2_log2"]
  dup1_median_NT_log2    <- median(dup1_NT_vec_log2)
  dup2_median_NT_log2    <- median(dup2_NT_vec_log2)

  x[, "GBA_rep1_log2_diff"]       <- x[, "GBA_rep1_log2"] - dup1_median_NT_log2
  x[, "GBA_rep2_log2_diff"]       <- x[, "GBA_rep2_log2"] - dup2_median_NT_log2

  # Normalize luminescent values
  lum_NT_vec                     <- x[are_NT, "Luminescence"]
  lum_norm_vec                   <- x[, "Luminescence"] / median(lum_NT_vec)
  x[, "Luminescence_normalized"] <- lum_norm_vec

  # Standardized Fluorescence by Luminescence
  x[, "GBA_rep1_Glo_standardized"] <- x[, "GBA_rep1_absolute"] /
                                      x[, "Luminescence_normalized"]

  x[, "GBA_rep2_Glo_standardized"] <- x[, "GBA_rep2_absolute"] /
                                      x[, "Luminescence_normalized"]

  # Log transform fluorescent values standardized by luminescence
  x[, "GBA_rep1_Glo_stand_log"] <- log2(x[, ("GBA_rep1_Glo_standardized")])
  x[, "GBA_rep2_Glo_stand_log"] <- log2(x[, ("GBA_rep2_Glo_standardized")])


  # Normalize Glo standardized values
  dup1_median_NT <- median(x[are_NT, "GBA_rep1_Glo_standardized"])
  dup2_median_NT <- median(x[are_NT, "GBA_rep2_Glo_standardized"])
  x[, "GBA_rep1_norm_Glo"] <- x[, "GBA_rep1_Glo_standardized"] / dup1_median_NT
  x[, "GBA_rep2_norm_Glo"] <- x[, "GBA_rep2_Glo_standardized"] / dup2_median_NT
  x[, "GBA_rep1_diff_Glo"] <- x[, "GBA_rep1_Glo_standardized"] - dup1_median_NT
  x[, "GBA_rep2_diff_Glo"] <- x[, "GBA_rep2_Glo_standardized"] - dup2_median_NT

  return(x)
})

GBA_df <- do.call(rbind.data.frame, c(split_df_list, make.row.names = FALSE))



# Calculate log-fold change -----------------------------------------------

GBA_df[, "GBA_rep1_logFC"]     <- log2(GBA_df[, "GBA_rep1_normalized"])
GBA_df[, "GBA_rep2_logFC"]     <- log2(GBA_df[, "GBA_rep2_normalized"])
GBA_df[, "GBA_rep1_Glo_logFC"] <- log2(GBA_df[, "GBA_rep1_norm_Glo"])
GBA_df[, "GBA_rep2_Glo_logFC"] <- log2(GBA_df[, "GBA_rep2_norm_Glo"])


# Calculate SSMD ----------------------------------------------------------

GBA_df[, "SSMD_MM_paired"]     <- Calculate_SSMD_var(GBA_df, "GBA_rep1_absolute")
GBA_df[, "SSMD_MM_paired_Glo"] <- Calculate_SSMD_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "SSMD_log2"]          <- Calculate_SSMD_var(GBA_df, "GBA_rep1_absolute", log2FC = TRUE)
GBA_df[, "SSMD_log2_Glo"]      <- Calculate_SSMD_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)

GBA_df[, "SSMD_test"]           <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_absolute")
GBA_df[, "SSMD_test_Glo"]       <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "SSMD_log2_test"]      <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_absolute", log2FC = TRUE)
GBA_df[, "SSMD_log2_Glo_test"]  <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)
GBA_df[, "SSMD_test_log2_diff"] <- Calculate_SSMD_var2(GBA_df, "GBA_rep1_log2_diff")

# Calculate T Score ----------------------------------------------------------

GBA_df[, "T_value"]         <- Calculate_T_var(GBA_df, "GBA_rep1_absolute")
GBA_df[, "T_value_Glo"]     <- Calculate_T_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "T_value_log"]     <- Calculate_T_var(GBA_df, "GBA_rep1_absolute", log2FC = TRUE)
GBA_df[, "T_value_Glo_log"] <- Calculate_T_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)


# Calculate P value ----------------------------------------------------------

GBA_df[, "P_value"]         <- Calculate_P_var(GBA_df, "GBA_rep1_absolute")
GBA_df[, "P_value_Glo"]     <- Calculate_P_var(GBA_df, "GBA_rep1_Glo_standardized")
GBA_df[, "P_value_log"]     <- Calculate_P_var(GBA_df, "GBA_rep1_absolute", log2FC = TRUE)
GBA_df[, "P_value_Glo_log"] <- Calculate_P_var(GBA_df, "GBA_rep1_Glo_standardized", log2FC = TRUE)



# Calculate hit strength --------------------------------------------------

GBA_df[, "Hit_strength"]     <- abs(rowMeans(log2(GBA_df[, c("GBA_rep1_normalized", "GBA_rep2_normalized")]))) *
                                -log10(GBA_df[, "P_value_log"])
GBA_df[, "Hit_strength_Glo"] <- abs(rowMeans(log2((GBA_df[, c("GBA_rep1_norm_Glo", "GBA_rep2_norm_Glo")])))) *
                                -log10(GBA_df[, "P_value_Glo_log"])




# Export Data -------------------------------------------------------------

write.csv(GBA_df,
          file = file.path(output_dir, "Tables", "GBA_complete.csv"),
          row.names = FALSE, quote = FALSE
          )



# Save data ---------------------------------------------------------------

save(GBA_df, file = file.path(r_data_dir, "03_analyse_data.RData"))
