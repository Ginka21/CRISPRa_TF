# 2021-12-27


# Load packages and source code -------------------------------------------

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "R_functions")
source(file.path(functions_dir, "01_calculating_scores.R"))



# Define folder path ------------------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
input_dir   <- file.path(project_dir, "2_input")
r_data_dir  <- file.path(project_dir, "3_R_objects")
output_dir  <- file.path(project_dir,"4_output")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_integrate_data.RData"))



# Normalize by non-targeting controls -------------------------------------

GBA_df[, "CellTiterGlo_foldNT"] <- NormPlates(GBA_df, "CellTiterGlo_raw", foldNT = TRUE)

for (i in 1:2) {
  ri <- paste0("_rep", i)

  GBA_df[, paste0("DeltaNT", ri)]             <- NormPlates(GBA_df, paste0("Raw", ri))
  GBA_df[, paste0("FoldNT", ri)]              <- NormPlates(GBA_df, paste0("Raw", ri), foldNT = TRUE)
  GBA_df[, paste0("PercActivation", ri)]      <- NormPlates(GBA_df, paste0("Raw", ri), percent_activation = TRUE)
  GBA_df[, paste0("Raw_log2", ri)]            <- log2(GBA_df[, paste0("Raw", ri)])
  GBA_df[, paste0("Log2FC", ri)]              <- NormPlates(GBA_df, paste0("Raw", ri), take_log2 = TRUE)
  GBA_df[, paste0("PercActivation_log2", ri)] <- NormPlates(GBA_df, paste0("Raw", ri), percent_activation = TRUE, take_log2 = TRUE)

  GBA_df[, paste0("Raw_Glo", ri)]                 <- GBA_df[, paste0("Raw", ri)] / GBA_df[, "CellTiterGlo_foldNT"]
  GBA_df[, paste0("DeltaNT_Glo", ri)]             <- NormPlates(GBA_df, paste0("Raw_Glo", ri))
  GBA_df[, paste0("FoldNT_Glo", ri)]              <- NormPlates(GBA_df, paste0("Raw_Glo", ri), foldNT = TRUE)
  GBA_df[, paste0("PercActivation_Glo", ri)]      <- NormPlates(GBA_df, paste0("Raw_Glo", ri), percent_activation = TRUE)
  GBA_df[, paste0("Raw_log2_Glo", ri)]            <- log2(GBA_df[, paste0("Raw_Glo", ri)])
  GBA_df[, paste0("Log2FC_Glo", ri)]              <- NormPlates(GBA_df, paste0("Raw_Glo", ri), take_log2 = TRUE)
  GBA_df[, paste0("PercActivation_log2_Glo", ri)] <- NormPlates(GBA_df, paste0("Raw_Glo", ri), percent_activation = TRUE, take_log2 = TRUE)
}

stripped_columns <- sub("_rep[12]", "", names(GBA_df))
GBA_df <- GBA_df[, order(match(stripped_columns, stripped_columns))]



# Calculate SSMD ----------------------------------------------------------

GBA_df[, "SSMD_deltaNT"]      <- Calculate_SSMD_var(GBA_df, "Raw_rep1")
GBA_df[, "SSMD_act"]          <- Calculate_SSMD_var(GBA_df, "Raw_rep1", percent_activation = TRUE)
GBA_df[, "SSMD_log2"]         <- Calculate_SSMD_var(GBA_df, "Raw_rep1", take_log2 = TRUE)
GBA_df[, "SSMD_act_log2"]     <- Calculate_SSMD_var(GBA_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE)

GBA_df[, "SSMD_deltaNT_Glo"]  <- Calculate_SSMD_var(GBA_df, "Raw_Glo_rep1")
GBA_df[, "SSMD_act_Glo"]      <- Calculate_SSMD_var(GBA_df, "Raw_Glo_rep1", percent_activation = TRUE)
GBA_df[, "SSMD_log2_Glo"]     <- Calculate_SSMD_var(GBA_df, "Raw_Glo_rep1", take_log2 = TRUE)
GBA_df[, "SSMD_act_log2_Glo"] <- Calculate_SSMD_var(GBA_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE)



# Calculate P value ----------------------------------------------------------

GBA_df[, "p_value_deltaNT"]      <- Calculate_P_var(GBA_df, "Raw_rep1")
GBA_df[, "p_value_act"]          <- Calculate_P_var(GBA_df, "Raw_rep1", percent_activation = TRUE)
GBA_df[, "p_value_log2"]         <- Calculate_P_var(GBA_df, "Raw_rep1", take_log2 = TRUE)
GBA_df[, "p_value_act_log2"]     <- Calculate_P_var(GBA_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE)

GBA_df[, "p_value_deltaNT_Glo"]  <- Calculate_P_var(GBA_df, "Raw_Glo_rep1")
GBA_df[, "p_value_act_Glo"]      <- Calculate_P_var(GBA_df, "Raw_Glo_rep1", percent_activation = TRUE)
GBA_df[, "p_value_log2_Glo"]     <- Calculate_P_var(GBA_df, "Raw_Glo_rep1", take_log2 = TRUE)
GBA_df[, "p_value_act_log2_Glo"] <- Calculate_P_var(GBA_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE)




# Calculate hit strength --------------------------------------------------

meanFC <- rowMeans(GBA_df[, c("Log2FC_rep1", "Log2FC_rep2")])
GBA_df[, "Hit_strength_deltaNT"]      <- meanFC * -log10(GBA_df[, "p_value_deltaNT"])
GBA_df[, "Hit_strength_act"]          <- meanFC * -log10(GBA_df[, "p_value_act"])
GBA_df[, "Hit_strength_log2"]         <- meanFC * -log10(GBA_df[, "p_value_log2"])
GBA_df[, "Hit_strength_act_log2"]     <- meanFC * -log10(GBA_df[, "p_value_act_log2"])

meanFC_Glo <- rowMeans(GBA_df[, c("Log2FC_Glo_rep1", "Log2FC_Glo_rep2")])
GBA_df[, "Hit_strength_deltaNT_Glo"]  <- meanFC_Glo * -log10(GBA_df[, "p_value_deltaNT_Glo"])
GBA_df[, "Hit_strength_act_Glo"]      <- meanFC_Glo * -log10(GBA_df[, "p_value_act_Glo"])
GBA_df[, "Hit_strength_log2_Glo"]     <- meanFC_Glo * -log10(GBA_df[, "p_value_log2_Glo"])
GBA_df[, "Hit_strength_act_log2_Glo"] <- meanFC_Glo * -log10(GBA_df[, "p_value_act_log2_Glo"])




# Export data -------------------------------------------------------------

write.csv(GBA_df,
          file = file.path(output_dir, "Tables", "GBA_complete.csv"),
          row.names = FALSE, quote = FALSE
          )



# Save data ---------------------------------------------------------------

save(GBA_df, file = file.path(r_data_dir, "03_analyse_data.RData"))


