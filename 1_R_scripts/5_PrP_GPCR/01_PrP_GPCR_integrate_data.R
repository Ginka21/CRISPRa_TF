# 2022-09-19



# Import packages and source code -----------------------------------------

library("readxl")



# Define folder path ------------------------------------------------------

project_dir       <- "~/R_projects/CRISPRa_TF"
input_dir         <- file.path(project_dir, "2_input")
rdata_dir         <- file.path(project_dir, "3_R_objects")
general_rdata_dir <- file.path(rdata_dir, "1_General")
PrP_rdata_dir     <- file.path(rdata_dir, "3_PrP", "GPCRa")
raw_data_dir      <- file.path(input_dir, "PrP_data")
layout_path       <- file.path(input_dir, "PrP_data", "LAYOUT", "GPCRa", "Layout.xlsx")



# Load data ---------------------------------------------------------------

load(file.path(input_dir, "CRISPRa_4sg_df.RData"))



# Read in data ------------------------------------------------------------

sheets_vec <- excel_sheets(layout_path)
layout_mat_list <- sapply(sheets_vec[1:4], function(x) {
  as.matrix(read_excel(layout_path, sheet = x)[1:16, 2:25])
}, simplify = FALSE)




# Define functions --------------------------------------------------------

OrderFileNames <- function(short_file_names) {
  file_splits <- strsplit(short_file_names, "_", fixed = TRUE)
  roman_numerals <- sapply(file_splits, "[[", 2)
  new_order <- order(as.integer(as.roman(roman_numerals)))
  return(new_order)
}

ReadInPrPData <- function(sub_folder, skip_lines) {
  folder_path <- path.expand(file.path(raw_data_dir, sub_folder))
  file_names_vec <- list.files(folder_path)
  file_names_vec <- grep("\\.csv$", file_names_vec, value = TRUE)
  paths_vec <- file.path(folder_path, file_names_vec)
  df_list <- lapply(paths_vec, function(x) {
    results_df <- read.csv(x, skip = skip_lines, stringsAsFactors = FALSE,
                           blank.lines.skip = FALSE)
    if (ncol(results_df) == 1) {
      results_df <- read.csv(x, skip = skip_lines, stringsAsFactors = FALSE,
                             blank.lines.skip = FALSE,
                             sep = ";"
                             )
    }
    return(results_df)
  })
  names(df_list) <- file_names_vec
  df_list <- df_list[OrderFileNames(names(df_list))]
  return(df_list)
}



# Read in data ------------------------------------------------------------

FRET_df_list <- ReadInPrPData("TR-FRET/GPCRa", skip_lines = 38)



# Process FRET values -----------------------------------------------------

ch1_mat_list <- lapply(FRET_df_list, function(x) {
  results_mat <- as.matrix(x[1:16, 2:25])
  mode(results_mat) <- "integer"
  colnames(results_mat) <- NULL
  return(results_mat)
})
ch2_mat_list <- lapply(FRET_df_list, function(x) {
  results_mat <- as.matrix(x[29:44, 2:25])
  mode(results_mat) <- "integer"
  colnames(results_mat) <- NULL
  return(results_mat)
})

odd_rows  <- seq(1, 15, by = 2)
even_rows <- seq(2, 16, by = 2)

FRET_P_vec <- vapply(seq_along(FRET_df_list), function(x) {
  (mean(ch1_mat_list[[x]][even_rows, 1]) - mean(ch1_mat_list[[x]][even_rows, 24])) /
  (mean(ch2_mat_list[[x]][even_rows, 1]) - mean(ch2_mat_list[[x]][even_rows, 24]))
}, numeric(1))

FRET_mat_list <- lapply(seq_along(FRET_df_list), function(x) {
  ch1_mat <- ch1_mat_list[[x]]
  ch2_mat <- ch2_mat_list[[x]]
  ch1_mat <- ch1_mat - mean(ch1_mat[odd_rows, 1])
  ch2_mat <- ch2_mat - mean(ch2_mat[even_rows, 24])
  fret_mat <- ch1_mat - (FRET_P_vec[[x]] * ch2_mat)
  fret_mat <- fret_mat - mean(fret_mat[odd_rows, 24])
  return(fret_mat)
})



# "Flatten" values from matrices to vectors -------------------------------

PrP_signal_vec_list <- lapply(FRET_mat_list, function(x) as.vector(t(x)))




# Prepare layout_df -------------------------------------------------------

layout_df_list <- lapply(seq_along(layout_mat_list), function(x) {
  data.frame(
    "Plate_number_384" = sub("Plate_", "", names(layout_mat_list)[[x]], fixed = TRUE),
    "Well_number_384"  = 1:384,
    "Plasmid_ID"       = as.vector(t(layout_mat_list[[x]])),
    stringsAsFactors   = FALSE
  )
})
layout_df <- do.call(rbind.data.frame,
                     c(layout_df_list,
                       stringsAsFactors = FALSE,
                       make.row.names = FALSE
                       ))


CRISPRa_4sg_df[, "Plasmid_name"] <- paste0(CRISPRa_4sg_df[, "Gene_symbol"],
                                           ifelse(is.na(CRISPRa_4sg_df[, "TSS_ID"]),
                                                  "",
                                                  paste0("_", CRISPRa_4sg_df[, "TSS_ID"])
                                                  )
                                           )

matches_vec <- match(layout_df[, "Plasmid_ID"], CRISPRa_4sg_df[, "Plasmid_name"])

use_columns <- c("Gene_symbol", "Entrez_ID", "TSS_ID", "Is_main_TSS")

layout_df <- data.frame(
  layout_df,
  CRISPRa_4sg_df[matches_vec, use_columns],
  "Is_NT_ctrl"  = layout_df[, "Plasmid_ID"] == "NT",
  "Is_pos_ctrl" = layout_df[, "Plasmid_ID"] == "PosCtrl",
  stringsAsFactors = FALSE,
  row.names = NULL
)




# Integrate measurement data with 'layout_df' -----------------------------

num_plates <- length(FRET_mat_list)
dup1_indices <- seq(from = 1, to = num_plates, by = 2)
dup2_indices <- seq(from = 2, to = num_plates, by = 2)
PrP_dup1_vec <- unlist(PrP_signal_vec_list[dup1_indices], use.names = FALSE)
PrP_dup2_vec <- unlist(PrP_signal_vec_list[dup2_indices], use.names = FALSE)

PrP_df <- data.frame(
  layout_df,
  "Raw_rep1" = PrP_dup1_vec,
  "Raw_rep2" = PrP_dup2_vec
)



# Save data ---------------------------------------------------------------

save(PrP_df, file = file.path(PrP_rdata_dir, "01_integrate_data.RData"))



