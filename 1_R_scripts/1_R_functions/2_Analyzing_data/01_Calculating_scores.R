# 2021-12-28



# Functions for data normalization ----------------------------------------

ObtainNTFactors <- function(input_df, use_column) {
  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  numeric_vec <- input_df[, use_column]
  mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
  are_own_NT <- input_df[, "Is_NT_ctrl"] & (input_df[, "Well_number_384"] %in% mat_384[, c(2, 22)])
  are_Tubingen_NT <- input_df[, "Is_NT_ctrl"] & !(are_own_NT)
  both_NT_plates <- unique(plate_numbers_vec[are_Tubingen_NT])
  median_factors <- vapply(both_NT_plates, function(x) {
    are_this_plate <- plate_numbers_vec == x
    median_own_NT <- median(numeric_vec[are_this_plate & are_own_NT])
    median_Tubingen_NT <- median(numeric_vec[are_this_plate & are_Tubingen_NT])
    NT_factor <- median_own_NT / median_Tubingen_NT
    return(NT_factor)
  }, numeric(1))
  names(median_factors) <- both_NT_plates
  return(median_factors)
}


ObtainAllNTFactors <- function(input_df, Glo = FALSE) {
  if (Glo) {
    median_factors_mat <- as.matrix(ObtainNTFactors(input_df, "CellTiterGlo_raw"))
  } else {
    median_factors_mat <- cbind(ObtainNTFactors(input_df, "Raw_rep1"),
                                ObtainNTFactors(input_df, "Raw_rep2")
                                )
    colnames(median_factors_mat) <- c("Raw_rep1", "Raw_rep2")
  }
  return(median_factors_mat)
}


NormPlates <- function(input_df,
                       use_column,
                       take_log2 = FALSE,
                       foldNT = FALSE,
                       percent_activation = FALSE,
                       correct_NT = FALSE
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

  if (correct_NT) {
    NT_factors <- ObtainAllNTFactors(input_df, Glo = use_column == "CellTiterGlo_raw")
    use_NT_factor <- median(NT_factors)
  }

  results_vec_list <- tapply(seq_along(numeric_vec), plate_numbers_vec, function(x) {

    sub_vec <- numeric_vec[x]
    sub_are_NT <- are_NT[x]
    sub_are_pos <- are_pos[x]

    median_NT <- median(sub_vec[sub_are_NT])
    if (correct_NT) {
      median_NT <- median_NT / use_NT_factor
    }

    if (foldNT && !(take_log2)) {
      sub_results <- sub_vec / median_NT
    } else {
      sub_results <- sub_vec - median_NT
      if (percent_activation) {
        sub_results <- sub_results / (median(sub_results[sub_are_pos]) - median(sub_results[sub_are_NT]))
      }
    }
    return(sub_results)
  })

  return(unlist(results_vec_list, use.names = FALSE))
}



# Functions for calculating test statistics for individual genes ----------

Calculate_SSMD <- function(input_df,
                           rep1_column,
                           t_score = FALSE,
                           plate_wise_NT_variance = TRUE,
                           ...
                           ) {

  rep2_column <- sub("_rep1", "_rep2", rep1_column, fixed = TRUE)

  norm_rep1 <- NormPlates(input_df, rep1_column, ...)
  norm_rep2 <- NormPlates(input_df, rep2_column, ...)

  are_NT <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  if (!(plate_wise_NT_variance)) {
    var_vec <- mapply(function(x, y) var(c(x, y)), norm_rep1[are_NT], norm_rep2[are_NT])
  }

  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))

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



# Functions for plate-level quality control -------------------------------

Calculate_Z_Prime <- function(sub_df, use_column, filter_NT = FALSE) {

  are_NT <- sub_df[, "Is_NT_ctrl"]
  if (filter_NT) {
    mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
    are_NT <- are_NT & (sub_df[, "Well_number_384"] %in% mat_384[, c(2, 22)])
  }
  are_pos <- sub_df[, "Is_pos_ctrl"]

  # Calculate Z' factor
  ## using means and standard deviation of controls
  NT_vec       <- sub_df[are_NT, use_column]
  posctrl_vec  <- sub_df[are_pos, use_column]
  mean_NT      <- mean(NT_vec)
  mean_posctrl <- mean(posctrl_vec)
  sd_NT        <- sd(NT_vec)
  sd_posctrl   <- sd(posctrl_vec)

  z_prime      <- (1 - (3 * (sd_NT + sd_posctrl)) / (mean_posctrl - mean_NT))
  return(z_prime)
}


Calculate_SSMD_ctrls <- function(sub_df, use_column, filter_NT = FALSE) {

  are_NT <- sub_df[, "Is_NT_ctrl"]
  if (filter_NT) {
    mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
    are_NT <- are_NT & (sub_df[, "Well_number_384"] %in% mat_384[, c(2, 22)])
  }
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


