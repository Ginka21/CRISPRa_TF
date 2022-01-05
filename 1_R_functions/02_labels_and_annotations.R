# 2121-12-29



# Define labels -----------------------------------------------------------

column_labels <- c(
  "GBA_rep1_absolute"         = "GBA activity (absolute values)",
  "GBA_rep1_diff"             = "GBA activity (difference from NT controls)",
  "GBA_rep1_normalized"       = "GBA activity (normalized to NT controls)",
  "GBA_rep1_log2"             = "GBA activity (log2 of absolute values)",
  "GBA_rep1_log2_diff"        = "GBA activity (log median NT substracted from log2 values)",
  "GBA_rep1_Glo_standardized" = "GBA activity (standardized using CellTiterGlo)",
  "GBA_rep1_Glo_stand_log"    = "GBA activity (log2 of standardized using CellTiterGlo)",

  "GBA_rep1_logFC"            = "GBA activity (log2 fold change)",
  "GBA_rep1_Glo_logFC"        = "GBA activity (log2 fold change of Glo standardized values",

  "Luminescence"              = "Cell viability (CellTiterGlo values)",
  "Luminescence_normalized"   = "Cell viability (normalized to NT controls)",

  "SSMD_MM_paired"            = "SSMD",
  "SSMD_MM_paired_Glo"        = "SSMD (standardized using CellTiterGlo)",
  "SSMD_log2"                 = "SSMD (using log-fold values)",
  "SSMD_log2_Glo"             = "SSMD (log-fold, standardized using CellTiterGlo)",

  "P_value"                   = "P values",
  "P_value_Glo"               = "P values (standardized using CellTiterGlo)",
  "P_value_log"               = "P values (using log-fold values)",
  "P_value_Glo_log"           = "P values (log-fold, standardized using CellTiterGlo)",

  "Hit_strength"              = "Hit strength",
  "Hit_strength_Glo"          = "Hit strength (standardized using CellTiterGlo)"


)