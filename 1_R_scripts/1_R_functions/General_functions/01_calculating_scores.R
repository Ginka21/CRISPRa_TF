# 2021-12-28


# Define functions --------------------------------------------------------

Calculate_Z_Prime <- function(sub_df, use_column) {

  are_NT <- sub_df[, "Is_NT_ctrl"]
  are_pos <- sub_df[, "Is_pos_ctrl"]

  # Calculate Z' factor
  ## Means and Standard Deviation of controls
  NT_vec       <- sub_df[are_NT, use_column]
  posctrl_vec  <- sub_df[are_pos, use_column]
  mean_NT      <- mean(NT_vec)
  mean_posctrl <- mean(posctrl_vec)
  sd_NT        <- sd(NT_vec)
  sd_posctrl   <- sd(posctrl_vec)

  z_prime      <- (1 - (3 * (sd_NT + sd_posctrl)) / (mean_posctrl - mean_NT))
  return(z_prime)
}


Calculate_SSMD_ctrls <- function(sub_df, use_column) {

  are_NT <- sub_df[, "Is_NT_ctrl"]
  are_pos <- sub_df[, "Is_pos_ctrl"]

  ## Means and variance of controls
  NT_vec       <- sub_df[are_NT, use_column]
  posctrl_vec  <- sub_df[are_pos, use_column]
  mean_NT      <- mean(NT_vec)
  mean_posctrl <- mean(posctrl_vec)
  var_NT       <- var(NT_vec)
  var_posctrl  <- var(posctrl_vec)

  SSMD_ctrl    <- (mean_posctrl - mean_NT) / (sqrt(var_posctrl + var_NT))
  return(SSMD_ctrl)
}



NormPlates <- function(input_df,
                       use_column,
                       take_log2 = FALSE,
                       foldNT = FALSE,
                       percent_activation = FALSE
                       ) {

  if (foldNT && percent_activation) {
    stop("The 'foldNT' and 'percent_activation' arguments are mutually exclusive!")
  }

  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  are_NT <- input_df[, "Is_NT_ctrl"]
  are_pos <- input_df[, "Is_pos_ctrl"]
  numeric_vec <- input_df[, use_column]
  if (take_log2) {
    numeric_vec <- log2(numeric_vec)
  }

  results_vec_list <- tapply(seq_along(numeric_vec), plate_numbers_vec, function(x) {

    sub_vec <- numeric_vec[x]
    sub_are_NT <- are_NT[x]
    sub_are_pos <- are_pos[x]

    median_NT <- median(sub_vec[sub_are_NT])

    if (foldNT && !(take_log2)) {
      sub_results <- sub_vec / median_NT
    } else {
      sub_results <- sub_vec - median_NT
      if (percent_activation) {
        sub_results <- sub_results / (median(sub_vec[sub_are_pos]) - median_NT)
      }
    }
    return(sub_results)
  })

  return(unlist(results_vec_list, use.names = FALSE))
}




Calculate_SSMD <- function(input_df,
                           rep1_column,
                           t_score = FALSE,
                           plate_wise_NT_variance = TRUE,
                           ...
                           ) {

  rep2_column <- sub("_rep1", "_rep2", rep1_column, fixed = TRUE)

  norm_rep1 <- NormPlates(input_df, rep1_column, ...)
  norm_rep2 <- NormPlates(input_df, rep2_column, ...)

  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  split_df_list <- split(input_df, plate_numbers_vec)

  are_NT <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")

  if (!(plate_wise_NT_variance)) {
    var_vec <- mapply(function(x, y) var(c(x, y)), norm_rep1[are_NT], norm_rep2[are_NT])
  }

  results_vec_list <- tapply(seq_along(are_NT), plate_numbers_vec, function(x) {

    sub_are_NT <- are_NT[x]
    rep1_vec <- norm_rep1[x]
    rep2_vec <- norm_rep2[x]

    if (plate_wise_NT_variance) {
      var_vec <- mapply(function(x, y) var(c(x, y)), rep1_vec[sub_are_NT], rep2_vec[sub_are_NT])
    }
    median_NT_var <- median(var_vec)

    results_vec <- vapply(seq_along(x), function(y) {
      delta_vec <- c(rep1_vec[[y]], rep2_vec[[y]])
      mean_diff <- mean(delta_vec)
      var_diff <- var(delta_vec)
      divisor <- (0.5 * var_diff) + (0.5 * median_NT_var)
      if (t_score) {
        divisor <- divisor / 2
      }
      return(mean_diff / sqrt(divisor))
    }, numeric(1))

  })

  return(unlist(results_vec_list, use.names = FALSE))
}


Calculate_T <- function(input_df, rep1_column, ...) {
  Calculate_SSMD(input_df, rep1_column, t_score = TRUE, ...)
}


Calculate_P <- function(input_df, rep1_column, ...) {
  t_values_vec <- Calculate_SSMD(input_df, rep1_column, ...)
  p_values_vec <- (2 * pt(abs(t_values_vec), 1, lower.tail = FALSE))
  return(p_values_vec)
}


