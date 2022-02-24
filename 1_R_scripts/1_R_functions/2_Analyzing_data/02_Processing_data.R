# 2022-01-17



# Define functions --------------------------------------------------------

NormalizeWithNTControls <- function(use_df, correct_NT = FALSE) {

  use_df[, "CellTiterGlo_foldNT"] <- NormPlates(use_df, "CellTiterGlo_raw", foldNT = TRUE, correct_NT = correct_NT)

  are_after <- seq_len(ncol(use_df)) > which(names(use_df) == "CellTiterGlo_raw")
  use_columns <- unique(c(names(use_df)[!(are_after)], "CellTiterGlo_foldNT",
                          names(use_df)[are_after]
                          ))
  use_df <- use_df[, use_columns]

  for (ri in paste0("_rep", 1:2)) {
    use_df[, paste0("DeltaNT", ri)]             <- NormPlates(use_df, paste0("Raw", ri), correct_NT = correct_NT)
    use_df[, paste0("FoldNT", ri)]              <- NormPlates(use_df, paste0("Raw", ri), foldNT = TRUE, correct_NT = correct_NT)
    use_df[, paste0("PercActivation", ri)]      <- NormPlates(use_df, paste0("Raw", ri), percent_activation = TRUE, correct_NT = correct_NT)
    use_df[, paste0("Raw_log2", ri)]            <- log2(use_df[, paste0("Raw", ri)])
    use_df[, paste0("Log2FC", ri)]              <- NormPlates(use_df, paste0("Raw", ri), take_log2 = TRUE, correct_NT = correct_NT)
    use_df[, paste0("PercActivation_log2", ri)] <- NormPlates(use_df, paste0("Raw", ri), percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)

    use_df[, paste0("Raw_Glo", ri)]                 <- use_df[, paste0("Raw", ri)] / use_df[, "CellTiterGlo_foldNT"]
    use_df[, paste0("DeltaNT_Glo", ri)]             <- NormPlates(use_df, paste0("Raw_Glo", ri), correct_NT = correct_NT)
    use_df[, paste0("FoldNT_Glo", ri)]              <- NormPlates(use_df, paste0("Raw_Glo", ri), foldNT = TRUE, correct_NT = correct_NT)
    use_df[, paste0("PercActivation_Glo", ri)]      <- NormPlates(use_df, paste0("Raw_Glo", ri), percent_activation = TRUE, correct_NT = correct_NT)
    use_df[, paste0("Raw_log2_Glo", ri)]            <- log2(use_df[, paste0("Raw_Glo", ri)])
    use_df[, paste0("Log2FC_Glo", ri)]              <- NormPlates(use_df, paste0("Raw_Glo", ri), take_log2 = TRUE, correct_NT = correct_NT)
    use_df[, paste0("PercActivation_log2_Glo", ri)] <- NormPlates(use_df, paste0("Raw_Glo", ri), percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)
  }

  stripped_columns <- sub("_rep[12]", "", names(use_df))
  use_df <- use_df[, order(match(stripped_columns, stripped_columns))]

  return(use_df)
}



RunSSMDStats <- function(use_df, correct_NT = FALSE) {

  ## Calculate SSMD
  use_df[, "SSMD_deltaNT"]      <- Calculate_SSMD(use_df, "Raw_rep1", correct_NT = correct_NT)
  use_df[, "SSMD_act"]          <- Calculate_SSMD(use_df, "Raw_rep1", percent_activation = TRUE, correct_NT = correct_NT)
  use_df[, "SSMD_log2"]         <- Calculate_SSMD(use_df, "Raw_rep1", take_log2 = TRUE, correct_NT = correct_NT)
  use_df[, "SSMD_act_log2"]     <- Calculate_SSMD(use_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)

  use_df[, "SSMD_deltaNT_Glo"]  <- Calculate_SSMD(use_df, "Raw_Glo_rep1", correct_NT = correct_NT)
  use_df[, "SSMD_act_Glo"]      <- Calculate_SSMD(use_df, "Raw_Glo_rep1", percent_activation = TRUE, correct_NT = correct_NT)
  use_df[, "SSMD_log2_Glo"]     <- Calculate_SSMD(use_df, "Raw_Glo_rep1", take_log2 = TRUE, correct_NT = correct_NT)
  use_df[, "SSMD_act_log2_Glo"] <- Calculate_SSMD(use_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)


  ## Calculate p value
  use_df[, "p_value_deltaNT"]      <- Calculate_P(use_df, "Raw_rep1", correct_NT = correct_NT)
  use_df[, "p_value_act"]          <- Calculate_P(use_df, "Raw_rep1", percent_activation = TRUE, correct_NT = correct_NT)
  use_df[, "p_value_log2"]         <- Calculate_P(use_df, "Raw_rep1", take_log2 = TRUE, correct_NT = correct_NT)
  use_df[, "p_value_act_log2"]     <- Calculate_P(use_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)

  use_df[, "p_value_deltaNT_Glo"]  <- Calculate_P(use_df, "Raw_Glo_rep1", correct_NT = correct_NT)
  use_df[, "p_value_act_Glo"]      <- Calculate_P(use_df, "Raw_Glo_rep1", percent_activation = TRUE, correct_NT = correct_NT)
  use_df[, "p_value_log2_Glo"]     <- Calculate_P(use_df, "Raw_Glo_rep1", take_log2 = TRUE, correct_NT = correct_NT)
  use_df[, "p_value_act_log2_Glo"] <- Calculate_P(use_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE, correct_NT = correct_NT)


  ## Calculate hit strength
  meanFC <- rowMeans(use_df[, c("Log2FC_rep1", "Log2FC_rep2")])
  use_df[, "Hit_strength_deltaNT"]      <- meanFC * -log10(use_df[, "p_value_deltaNT"])
  use_df[, "Hit_strength_act"]          <- meanFC * -log10(use_df[, "p_value_act"])
  use_df[, "Hit_strength_log2"]         <- meanFC * -log10(use_df[, "p_value_log2"])
  use_df[, "Hit_strength_act_log2"]     <- meanFC * -log10(use_df[, "p_value_act_log2"])

  meanFC_Glo <- rowMeans(use_df[, c("Log2FC_Glo_rep1", "Log2FC_Glo_rep2")])
  use_df[, "Hit_strength_deltaNT_Glo"]  <- meanFC_Glo * -log10(use_df[, "p_value_deltaNT_Glo"])
  use_df[, "Hit_strength_act_Glo"]      <- meanFC_Glo * -log10(use_df[, "p_value_act_Glo"])
  use_df[, "Hit_strength_log2_Glo"]     <- meanFC_Glo * -log10(use_df[, "p_value_log2_Glo"])
  use_df[, "Hit_strength_act_log2_Glo"] <- meanFC_Glo * -log10(use_df[, "p_value_act_log2_Glo"])

  return(use_df)
}




