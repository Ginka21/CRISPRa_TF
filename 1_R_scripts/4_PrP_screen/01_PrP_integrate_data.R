# 2021-12-22



# Define folder path ------------------------------------------------------

project_dir       <- "~/R_projects/CRISPRa_TF"
input_dir         <- file.path(project_dir, "2_input")
rdata_dir         <- file.path(project_dir, "3_R_objects")
general_rdata_dir <- file.path(rdata_dir, "1_General")
PrP_rdata_dir     <- file.path(rdata_dir, "3_PrP")
raw_dir           <- file.path(input_dir, "PrP_data")



# Load data ---------------------------------------------------------------

load(file.path(general_rdata_dir, "01_convert_plate_layouts.RData"))



# Read in data ------------------------------------------------------------

FRET_folders <- list.files(file.path(raw_dir, "TR-FRET"), full.names = TRUE)
FRET_file_paths <- unlist(lapply(FRET_folders, function(x) {
  grep("\\.csv$", list.files(x, full.names = TRUE), value = TRUE)
}))
FRET_df_list <- lapply(FRET_file_paths, read.csv, skip = 38,
                       stringsAsFactors = FALSE, header = TRUE
                       )

Glo_folders <- list.files(file.path(raw_dir, "CELL TITER GLO"), full.names = TRUE)
Glo_file_paths <- unlist(lapply(Glo_folders, function(x) {
  grep("\\.csv$", list.files(x, full.names = TRUE), value = TRUE)
}))
Glo_df_list <- lapply(Glo_file_paths, read.csv, skip = 9, stringsAsFactors = FALSE)




# Define functions --------------------------------------------------------

TidyFileNames <- function(file_paths_vec) {
  file_splits <- strsplit(file_paths_vec, "5000_", fixed = TRUE)
  file_splits <- strsplit(sapply(file_splits, "[[", 2), "_PRP_", fixed = TRUE)
  results_vec <- sapply(file_splits, "[[", 1)
  return(results_vec)
}

OrderFileNames <- function(short_file_names) {
  file_splits <- strsplit(short_file_names, "_", fixed = TRUE)
  roman_numerals <- sapply(file_splits, "[[", 2)
  new_order <- order(as.integer(as.roman(roman_numerals)))
  return(new_order)
}



# Re-order data -----------------------------------------------------------

FRET_short_names <- TidyFileNames(FRET_file_paths)
names(FRET_df_list) <- FRET_short_names
FRET_df_list <- FRET_df_list[OrderFileNames(FRET_short_names)]

Glo_mat_list <- lapply(Glo_df_list, function(x) {
  results_mat <- as.matrix(x[1:16, 2:25])
  mode(results_mat) <- "integer"
  return(results_mat)
})
Glo_short_names <- TidyFileNames(Glo_file_paths)
names(Glo_mat_list) <- Glo_short_names
Glo_mat_list <- Glo_mat_list[OrderFileNames(Glo_short_names)]




# Process FRET values -----------------------------------------------------

ch1_mat_list <- lapply(FRET_df_list, function(x) {
  results_mat <- as.matrix(x[1:16, 2:25])
  mode(results_mat) <- "integer"
  colnames(results_mat) <- NULL
  return(results_mat)
})
ch2_mat_list <- lapply(FRET_df_list, function(x) {
  results_mat <- as.matrix(x[26:41, 2:25])
  mode(results_mat) <- "integer"
  colnames(results_mat) <- NULL
  return(results_mat)
})

every_2nd <- seq(2, 16, by = 2)
every_1st <- seq(1, 15, by = 2)

FRET_P_vec <- vapply(1:24, function(x) {
  (mean(ch1_mat_list[[x]][every_2nd, 1]) - mean(ch1_mat_list[[x]][every_2nd, 24])) /
  (mean(ch2_mat_list[[x]][every_2nd, 1]) - mean(ch2_mat_list[[x]][every_2nd, 24]))
}, numeric(1))

FRET_mat_list <- lapply(1:24, function(x) {
  ch1_mat <- ch1_mat_list[[x]]
  ch2_mat <- ch2_mat_list[[x]]
  ch1_mat <- ch1_mat - mean(ch1_mat[every_1st, 1])
  ch2_mat <- ch2_mat - mean(ch2_mat[every_2nd, 24])
  fret_mat <- ch1_mat - (FRET_P_vec[[x]] * ch2_mat)
  fret_mat <- fret_mat - mean(fret_mat[every_1st, 24])
  return(fret_mat)
})




# "Flatten" values from matrices to vectors -------------------------------

PrP_signal_vec_list <- lapply(FRET_mat_list, function(x) {
  as.vector(t(as.matrix(x)))
})

lum_signal_vec_list <- lapply(Glo_mat_list, function(x) {
  as.vector(t(x))
})



# Integrate signal values into layout_df ----------------------------------

dup1_indices <- seq(from = 1, to = 23, by = 2)
dup2_indices <- seq(from = 2, to = 24, by = 2)
PrP_dup1_vec <- unlist(PrP_signal_vec_list[dup1_indices], use.names = FALSE)
PrP_dup2_vec <- unlist(PrP_signal_vec_list[dup2_indices], use.names = FALSE)
lum_vec      <- unlist(lum_signal_vec_list, use.names = FALSE)
lum_vec      <- as.numeric(lum_vec)

PrP_df <- data.frame(
  layout_df,
  "CellTiterGlo_raw" = lum_vec,
  "Raw_rep1" = PrP_dup1_vec,
  "Raw_rep2" = PrP_dup2_vec
)


# Save data ---------------------------------------------------------------

save(PrP_df, file = file.path(PrP_rdata_dir, "01_integrate_data.RData"))



