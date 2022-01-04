# 2121-12-28

# Define functions --------------------------------------------------------

Calculate_SSMD <- function(numeric_vec) {
  mean_difference    <- mean(numeric_vec)
  standard_deviation <- sd(numeric_vec)
  SSMD               <- mean_difference / standard_deviation
  return(SSMD)
}

Calculate_T <- function(numeric_vec) {
  mean_difference <- mean(numeric_vec)
  standard_deviation <- sd(numeric_vec)
  t_score <- (mean_difference / (standard_deviation / sqrt(length(numeric_vec))))
  return(t_score)
}

Calculate_P <- function(numeric_vec) {
  t_value <- Calculate_T(numeric_vec)
  p_value <- (2 * pt(abs(t_value), 1, lower.tail = FALSE))
  return(p_value)
}


Calculate_Z_Prime <- function(sub_df, use_column) {
  are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  are_posctrl <- sub_df[, "Target_flag"] %in% c("Pos. control")

  # Calculate Z' factor
  ## Means and Standard Deviation of controls
  NT_vec <- sub_df[are_NT, use_column]
  posctrl_vec <- sub_df[are_posctrl, use_column]
  mean_NT <- mean(NT_vec)
  mean_posctrl <- mean(posctrl_vec)
  sd_NT <- sd(NT_vec)
  sd_posctrl <- sd(posctrl_vec)

  z_prime <- (1 - (3 * (sd_NT + sd_posctrl)) / (mean_posctrl - mean_NT))
  return(z_prime)
}



Calculate_SSMD_var <- function(input_df, rep1_column, t_score = FALSE, log2FC = FALSE) {

  rep2_column <- sub("_rep1", "_rep2", rep1_column, fixed = TRUE)

  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  split_df_list <- split(input_df, plate_numbers_vec)

  results_vec_list <- lapply(split_df_list, function(sub_df) {
    are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
    num_NT <- sum(are_NT)

    if (log2FC) {
      rep1_diff_vec <- log2(sub_df[, rep1_column] / median(sub_df[are_NT, rep1_column]))
      rep2_diff_vec <- log2(sub_df[, rep2_column] / median(sub_df[are_NT, rep2_column]))
    } else {
      rep1_diff_vec <- sub_df[, rep1_column] - median(sub_df[are_NT, rep1_column])
      rep2_diff_vec <- sub_df[, rep2_column] - median(sub_df[are_NT, rep2_column])
    }

    results_vec <- vapply(seq_len(nrow(sub_df)), function(x) {
      delta_vec <- c(rep1_diff_vec[[x]], rep2_diff_vec[[x]])
      mean_diff <- mean(delta_vec)
      var_diff <- var(delta_vec)
      var_diff_NT <- var(c(rep1_diff_vec[are_NT], rep2_diff_vec[are_NT]))
      if (t_score) {
        var_diff <- var_diff / 2
        var_diff_NT <- var_diff_NT / num_NT
      }
      return(mean_diff / sqrt(var_diff + var_diff_NT))
    }, numeric(1))

  })

  return(unlist(results_vec_list, use.names = FALSE))
}


Calculate_T_var <- function(input_df, rep1_column, log2FC = FALSE) {
  Calculate_SSMD_var(input_df, rep1_column, t_score = TRUE, log2FC = log2FC)
}

Calculate_P_var <- function(input_df, rep1_column, log2FC = FALSE) {
  t_values_vec <- Calculate_SSMD_var(input_df, rep1_column, t_score = TRUE, log2FC = log2FC)
  p_values_vec <- (2 * pt(abs(t_values_vec), 1, lower.tail = FALSE))
  return(p_values_vec)
}



