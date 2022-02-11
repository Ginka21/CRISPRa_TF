# 2021-12-22


# Define folder path ------------------------------------------------------

project_dir       <- "~/R_projects/CRISPRa_TF"
input_dir         <- file.path(project_dir, "2_input")
rdata_dir         <- file.path(project_dir, "3_R_objects")
general_rdata_dir <- file.path(rdata_dir, "1_General")
GBA_rdata_dir     <- file.path(rdata_dir, "2_GBA")
raw_data_dir      <- file.path(input_dir, "GBA_data")



# Define functions --------------------------------------------------------

SortByRoman <- function(char_vec) {
  char_splits <- strsplit(char_vec, "[-_]")
  roman_vec <- vapply(char_splits, function(x) {
    are_roman <- x %in% as.character(as.roman(1:12))
    return(x[are_roman])
  }, "")
  numeric_roman_vec <- as.integer(as.roman(roman_vec))
  results_vec <- char_vec[order(numeric_roman_vec)]
  return(results_vec)
}



# Load data ---------------------------------------------------------------

load(file.path(general_rdata_dir, "01_convert_plate_layouts.RData"))



# Read in data ------------------------------------------------------------

raw_files_vec <- list.files(raw_data_dir)
are_luminescence <- grepl("Luminescence", raw_files_vec, fixed = TRUE)
GBA_files_vec <- raw_files_vec[!are_luminescence]
GBA_files_vec <- SortByRoman(GBA_files_vec)

GBA_df_list <- sapply(GBA_files_vec, function(x) {
  read.csv(file.path(raw_data_dir, x), header = FALSE)
}, simplify = FALSE)

GBA_signal_vec_list <- lapply(GBA_df_list, function(x) {
  as.vector(t(as.matrix(x)))
})

luminescence_files_vec <- raw_files_vec[are_luminescence]
luminescence_files_vec <- SortByRoman(luminescence_files_vec)
lum_df_list <- sapply(luminescence_files_vec, function(x) {
  read.csv(file.path(raw_data_dir, x), skip = 9)[1:16, 2:25]
}, simplify = FALSE)

lum_signal_vec_list <- lapply(lum_df_list, function(x) {
  as.vector(t(as.matrix(x)))
})



# Exclude problematic control wells from 'layout_df' ----------------------
## These were known a priori (dispensing issues!)

columns_1to4 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)[, 1:4]
are_problematic <- (layout_df[, "Plate_number_384"] == "X") &
                   (layout_df[, "Well_number_384"] %in% columns_1to4)
layout_df[, "Is_problematic"] <- are_problematic
layout_df[are_problematic, "Is_NT_ctrl"] <- FALSE
layout_df[are_problematic, "Is_pos_ctrl"] <- FALSE



# Integrate measurement data with 'layout_df' -----------------------------

dup1_indices <- seq(from = 1, to = 23, by = 2)
dup2_indices <- seq(from = 2, to = 24, by = 2)
GBA_dup1_vec <- unlist(GBA_signal_vec_list[dup1_indices], use.names = FALSE)
GBA_dup2_vec <- unlist(GBA_signal_vec_list[dup2_indices], use.names = FALSE)
lum_vec      <- unlist(lum_signal_vec_list, use.names = FALSE)
lum_vec      <- as.numeric(lum_vec)

GBA_df <- data.frame(
  layout_df,
  "CellTiterGlo_raw" = lum_vec,
  "Raw_rep1" = GBA_dup1_vec,
  "Raw_rep2" = GBA_dup2_vec
)



# Save data ---------------------------------------------------------------

save(GBA_df, file = file.path(GBA_rdata_dir, "01_integrate_data.RData"))



