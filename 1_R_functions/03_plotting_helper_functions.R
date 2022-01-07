# 2022-01-03



# Define functions --------------------------------------------------------

ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}


VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("gh")), use_expression, expression(phantom("gh")))
  return(ConcatenateExpressions(my_list, my_sep = ""))
}

StripExpression <- function(my_expression) {
  if (is.character(my_expression)) {
    literal_string <- paste0("\"", capture.output(cat(my_expression)), "\"")
  } else {
    literal_string <- paste0(capture.output(my_expression), collapse = "")
    if (substr(literal_string, 1L, 11L) == "expression(") {
      literal_string <- substr(literal_string, 12L, nchar(literal_string) - 1L)
    }
  }
  return(literal_string)
}

Darken <- function(color, factor = 1.4) {
  # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}


palify_cache_101 <- list()
Palify <- function(myhex, fraction_pale = 0.5) {
  if (myhex %in% names(palify_cache_101)) {
    color_vec <- palify_cache_101[[myhex]]
  } else {
    color_vec <- colorRampPalette(c(myhex, "#FFFFFF"))(101)
    palify_cache_101[[myhex]] <- color_vec
    assign("palify_cache_101", palify_cache_101, envir = globalenv())
  }
  color_vec[[round(fraction_pale * 100) + 1]]
}


DrawSideLegend <- function(labels_list,
                           use_colors,
                           use_pch = 16,
                           use_point_size = 1.2
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing the legend
  y_mid <- 0.5
  small_gap <- diff(grconvertY(c(0, 1.25), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * 1.75

  if (all(lengths(labels_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(labels_list))
    are_first <- rep(TRUE, length(labels_list))
  } else {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  x_text  <- 1 + diff(grconvertX(c(0, 0.75), from = "lines", to = "npc"))
  x_point <- 1 + diff(grconvertX(c(0, 0.9), from = "lines", to = "npc"))

  ## Draw the legend
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))

  points(x   = rep(grconvertX(x = x_point, from = "npc", to = "user"), length(labels_list)),
         y   = tapply(y_pos, groups_vec, mean),
         cex = use_point_size,
         pch = use_pch,
         col = use_colors,
         xpd = NA
         )

  return(invisible(NULL))
}

