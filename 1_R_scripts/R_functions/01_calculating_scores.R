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
  NT_vec       <- sub_df[are_NT, use_column]
  posctrl_vec  <- sub_df[are_posctrl, use_column]
  mean_NT      <- mean(NT_vec)
  mean_posctrl <- mean(posctrl_vec)
  sd_NT        <- sd(NT_vec)
  sd_posctrl   <- sd(posctrl_vec)

  z_prime      <- (1 - (3 * (sd_NT + sd_posctrl)) / (mean_posctrl - mean_NT))
  return(z_prime)
}

Calculate_SSMD_ctrls <- function(sub_df, use_column) {
  are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  are_posctrl <- sub_df[, "Target_flag"] %in% c("Pos. control")


  ## Means and Variance of controls
  NT_vec       <- sub_df[are_NT, use_column]
  posctrl_vec  <- sub_df[are_posctrl, use_column]
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
  are_NT <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
  are_pos <- input_df[, "Target_flag"] %in% "Pos. control"
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







Calculate_SSMD_var <- function(input_df,
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





#
#
#
#
# Calculate_SSMD_var <- function(input_df,
#                                rep1_column,
#                                t_score = FALSE,
#                                log2FC = FALSE,
#                                percent_activation = FALSE
#                                ) {
#
#   rep2_column <- sub("_rep1", "_rep2", rep1_column, fixed = TRUE)
#
#   plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
#   split_df_list <- split(input_df, plate_numbers_vec)
#
#   results_vec_list <- lapply(split_df_list, function(sub_df) {
#     are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
#     rep1_vec <- sub_df[, rep1_column]
#     rep2_vec <- sub_df[, rep2_column]
#     if (log2FC) {
#       rep1_vec <- log2(rep1_vec)
#       rep2_vec <- log2(rep2_vec)
#     }
#     rep1_diff_vec <- rep1_vec - median(rep1_vec[are_NT])
#     rep2_diff_vec <- rep2_vec - median(rep2_vec[are_NT])
#
#     if (percent_activation) {
#       rep1_diff_vec <- rep1_diff_vec / (median(rep1_diff_vec[are_pos]) - median(rep1_vec[are_NT]))
#       rep2_diff_vec <- rep2_diff_vec / (median(rep2_diff_vec[are_pos]) - median(rep2_vec[are_NT]))
#     }
#
#     results_vec <- vapply(seq_len(nrow(sub_df)), function(x) {
#       delta_vec <- c(rep1_diff_vec[[x]], rep2_diff_vec[[x]])
#       mean_diff <- mean(delta_vec)
#       var_diff <- var(delta_vec)
#       divisor <- (0.5 * var_diff) + (0.5 * median_NT_var)
#       if (t_score) {
#         divisor <- divisor / 2
#       }
#       return(mean_diff / sqrt(divisor))
#     }, numeric(1))
#
#   })
#
#
#
#   are_NT <- input_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
#   var_vec <- mapply(function(x, y) var(c(x, y)),
#                     input_df[are_NT, rep1_column],
#                     input_df[are_NT, rep2_column]
#                     )
#   median_NT_var <- median(var_vec)
#
#   results_vec_list <- lapply(split_df_list, function(sub_df) {
#     are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
#     are_pos <- sub_df[, "Target_flag"] %in% "Pos. control"
#     num_NT <- sum(are_NT)
#
#     rep1_vec <- sub_df[, rep1_column]
#     rep2_vec <- sub_df[, rep2_column]
#
#     if (log2FC) {
#       rep1_vec <- log2(rep1_vec)
#       rep2_vec <- log2(rep2_vec)
#     }
#
#     rep1_diff_vec <- rep1_vec - median(rep1_vec[are_NT])
#     rep2_diff_vec <- rep2_vec - median(rep2_vec[are_NT])
#
#     if (percent_activation) {
#       rep1_diff_vec <- rep1_diff_vec / (median(rep1_diff_vec[are_pos]) - median(rep1_vec[are_NT]))
#       rep2_diff_vec <- rep2_diff_vec / (median(rep2_diff_vec[are_pos]) - median(rep2_vec[are_NT]))
#     }
#
#     results_vec <- vapply(seq_len(nrow(sub_df)), function(x) {
#       delta_vec <- c(rep1_diff_vec[[x]], rep2_diff_vec[[x]])
#       mean_diff <- mean(delta_vec)
#       var_diff <- var(delta_vec)
#       divisor <- (0.5 * var_diff) + (0.5 * median_NT_var)
#       if (t_score) {
#         divisor <- divisor / 2
#       }
#       return(mean_diff / sqrt(divisor))
#     }, numeric(1))
#
#   })
#
#   return(unlist(results_vec_list, use.names = FALSE))
# }
#
#




Calculate_T_var <- function(input_df, rep1_column, ...) {
  Calculate_SSMD_var(input_df, rep1_column, t_score = TRUE, ...)
}

Calculate_P_var <- function(input_df, rep1_column, ...) {
  t_values_vec <- Calculate_SSMD_var(input_df, rep1_column, ...)
  p_values_vec <- (2 * pt(abs(t_values_vec), 1, lower.tail = FALSE))
  return(p_values_vec)
}


# Calculate_SSMD_var2 <- function(input_df, rep1_column, log2FC = FALSE) {
#
#   rep2_column <- sub("_rep1", "_rep2", rep1_column, fixed = TRUE)
#
#   plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
#   split_df_list <- split(input_df, plate_numbers_vec)
#
#   results_vec_list <- lapply(split_df_list, function(sub_df) {
#     are_NT <- sub_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
#
#     if (log2FC) {
#       rep1_diff_vec <- log2(sub_df[, rep1_column] / median(sub_df[are_NT, rep1_column]))
#       rep2_diff_vec <- log2(sub_df[, rep2_column] / median(sub_df[are_NT, rep2_column]))
#       NT1_diff_vec  <- log2(sub_df[are_NT, rep1_column] / median(sub_df[are_NT, rep1_column]))
#       NT2_diff_vec  <- log2(sub_df[are_NT, rep2_column] / median(sub_df[are_NT, rep2_column]))
#     } else {
#       rep1_diff_vec <- sub_df[, rep1_column] - median(sub_df[are_NT, rep1_column])
#       rep2_diff_vec <- sub_df[, rep2_column] - median(sub_df[are_NT, rep2_column])
#       NT1_diff_vec  <- sub_df[are_NT, rep1_column] - median(sub_df[are_NT, rep1_column])
#       NT2_diff_vec  <- sub_df[are_NT, rep2_column] - median(sub_df[are_NT, rep2_column])
#     }
#
#       results_vec  <- vapply(seq_len(nrow(sub_df)), function(x) {
#       delta_vec    <- c(rep1_diff_vec[[x]], rep2_diff_vec[[x]])
#       mean_diff    <- mean(delta_vec)
#       delta_NT_vec <- c(rep1_diff_vec[are_NT], rep2_diff_vec[are_NT])
#       mean_NT_diff <- mean(delta_NT_vec)
#       var_diff     <- var(delta_vec)
#       var_diff_NT  <- var(c(rep1_diff_vec[are_NT], rep2_diff_vec[are_NT]))
#
#       return((mean_diff - mean_NT_diff)/ sqrt(var_diff + var_diff_NT))
#     }, numeric(1))
#
#   })
#
#   return(unlist(results_vec_list, use.names = FALSE))
# }


