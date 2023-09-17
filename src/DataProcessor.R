# @date: 2023-09-17

library(tools)
source("DataHandler.R")

deconvolute_multidrug <- function(
  drug_printer_file, output_dir
) {
  drug_printer_data <- read_drug_printer(drug_printer_file)
  data_dir <- dirname(drug_printer_file)
  data_files_list <- list.files(data_dir)
  plate_ids <- unique(drug_printer_data$plateID)
  viabilities <- list()
  for (plate in plate_ids) {
    plate_reader_file <- subset(data_files_list,
            grepl(plate, data_files_list) & !endsWith(data_files_list, ".xlsx"))
    if (length(plate_reader_file) > 1) {
      print(plate_reader_file)
      print("Multiple files found for a single plate. Check file names.")
      q("no")
    }
    plate_drug_printer_data <- subset(drug_printer_data, plateID == plate)
  }
}