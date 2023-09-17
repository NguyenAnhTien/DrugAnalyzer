# @date: 2023-09-17
# @desc:
#     1) This file contains functions for loading data from files

library(readxl)
source("Math.R")
source("DataFrame.R")
source("StringHandler.R")


read_drug_printer <- function(
  drug_printer_file
) {
  # @desc:
  #     1) Read the drug printer file and store into a DataFrame object
  #     2) Extract the relevant columns (skip irrelevant columns)
  #     3) Rename the column names
  #     4) Removing rows without a drug name
  #     5) Split the drug name by space
  #     6) Round the concentration to 4 significant figures
  #     7) Convert the plate ID to character
  # @params:
  #     1) drug_printer_file: The path to the drug printer file
  # @return:
  #     1) extracted_data: The DataFrame object containing the extracted data
  drug_printer_df <- read_excel(drug_printer_file, sheet = "Tabular detail")
  extracted_data <- drug_printer_df[c("Plate ID", "Dispensed\r\nwell",
                          "Dispensed\r\nrow", "Dispensed\r\ncol", "Fluid name",
                          "Dispensed concentration")]
  column_names_mapping <- c("plateID", "well", "row", "column", "drug",
                            "concentration")
  extracted_data <- rename_columns(extracted_data, column_names_mapping)
  extracted_data <- extracted_data[!is.na(extracted_data$drug), ]
  extracted_data$drug <- sapply(extracted_data$drug, split_string)
  extracted_data$concentration <- sapply(extracted_data$concentration,
                                         round_sig, sig = 4)
  extracted_data$plateID <- as.character(extracted_data$plateID)
  return(extracted_data)
}

read_plate <- function(
  plate_file
) {
  raw_data <- read_plate_file(plate_file)
  raw_data <- raw_data[13:28]
  rows <- LETTERS[1:16]
  columns <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "23", "24")
  process_data <- split_data(raw_data, rows, columns)
  plate_data_dict <- create_plate_data_dict(process_data, rows, columns)
  plate_data <- create_plate_data_frame(plate_data_dict)
  return(plate_data)
}

read_plate_file <- function(
  plate_file
) {
  file <- file(plate_file, "r")
  raw_data <- readLines(file, n = -1)
  close(file)
  return(raw_data)
}

split_data <- function(
  raw_data, row_names, col_names
) {
  rows <- unlist(strsplit(raw_data, "\n"))
  process_data <- lapply(rows, function(row) {
    as.integer(unlist(strsplit(row, ";"))[-1])
  })
  return(process_data)
}

create_plate_data_dict <- function(
  process_data, row_names, column_names
) {
  plate_data_dict <- list()
  for (i_idx in seq_along(row_names)) {
    for (j_idx in seq_along(column_names)) {
      reading_code <- paste0(row_names[i_idx], column_names[j_idx])
      reading <- process_data[[i_idx]][[j_idx]]
      plate_data_dict[[reading_code]] <- reading
    }
  }
  return(plate_data_dict)
}

create_plate_data_frame <- function(
  plate_data_dict
) {
  plate_data <- data.frame(well = names(plate_data_dict),
                           result = unlist(plate_data_dict))
  colnames(plate_data) <- c("well", "result")
  return(plate_data)
}

