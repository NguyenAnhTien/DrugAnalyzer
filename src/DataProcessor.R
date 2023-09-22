# @date: 2023-09-17

library(tools)
source("DataHandler.R")
source("DataMatcher.R")
source("DrugDataHandler.R")

deconvolute_multidrug <- function(
  drug_printer_file, output_dir
) {
  drug_printer_data <- read_drug_printer(drug_printer_file)
  data_dir <- dirname(drug_printer_file)
  data_files_list <- list.files(data_dir)
  plate_ids <- unique(drug_printer_data$plateID)
  viabilities <- list()
  for (plate in plate_ids) {
    plate_data_file <- subset(data_files_list,
            grepl(plate, data_files_list) & !endsWith(data_files_list, ".xlsx"))
    if (length(plate_data_file) > 1) {
      print(plate_data_file)
      print("Multiple files found for a single plate. Check file names.")
      q("no")
    }
    plate_drug_printer_data <- subset(drug_printer_data, plateID == plate)
    plate_data_file <- file.path(data_dir, plate_data_file)
    plate_data <- read_plate(plate_data_file)
    drug_combinations <- cal_drug_combinations(plate_drug_printer_data)
    combo_wells <- cal_drug_combo_wells(
                                        plate_drug_printer_data,
                                        drug_combinations)
    for (drug_combo in drug_combinations) {
      filtered_drug_printer_data <- filter_data_by_drug_combo(
                                                       plate_drug_printer_data,
                                                       combo_wells, drug_combo)
      matched_data <- match_data(filtered_drug_printer_data, plate_data)
      drug_data_table <- create_data_table(matched_data)
      drug_screen_table <- save_drug_data(plate, matched_data, drug_data_table,
                                                                    output_dir)
    }
  }
  return(drug_screen_table)
}

filter_data_by_drug_combo <- function(
  plate_drug_printer_data, combo_wells, drug_combo
) {
  combo_wells_for_combo <- combo_wells[[toString(drug_combo)]]
  filtered_data <- subset(plate_drug_printer_data,
                          well %in% combo_wells_for_combo)
  return(filtered_data)
}