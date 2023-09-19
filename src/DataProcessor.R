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
      result <- create_data_table(matched_data)
      return(result)
    }
  }
}

create_data_table <- function(
  matched_data
) {
  matched_wells <- matched_data[[1]]
  drug_conc_ranges <- matched_data[[2]]
  # Extract drug names and remove unwanted drugs
  drug_names <- names(drug_conc_ranges)
  drug_names <- drug_names[!(drug_names %in% c("DMSO", "a+Tw", "a+B"))]
  
  # Ensure we're working with only two drugs
  if (length(drug_names) > 2) {
    stop("Only two drug combinations can be processed.")
  }
  
  # Extract the two drug names
  drug_name1 <- drug_names[1]
  drug_name2 <- drug_names[2]
  
  # Initialize data frames for mean, n, and standard deviation
  col_index <- drug_conc_ranges[[drug_name1]]
  row_index <- drug_conc_ranges[[drug_name2]]
  drug_screen_table_dict <- list()
  
  # Iterate over each well and populate the data frames
  for (well in names(matched_wells)) {
    drug_info <- matched_wells[[well]]$drug_info
    ind <- 0
    col <- 0
    
    for (drug in names(drug_info)) {
      if (drug != "DMSO" && drug != "a+Tw" && drug != "a+B") {
        if (drug == drug_name1) {
          col <- as.numeric(drug_info[[drug]])
        } else if (drug == drug_name2) {
          ind <- as.numeric(drug_info[[drug]])
        }
      }
    }
    
    # Add the reading to the appropriate cell in the data frames
    if (!is.null(col) && !is.null(ind)) {
      key <- paste(ind, col, sep = "_")
      if (!is.null(matched_wells[[well]]$reading)) {
        if (is.null(drug_screen_table_dict[[key]])) {
          drug_screen_table_dict[[key]] <- c(matched_wells[[well]]$reading)
        } else {
          drug_screen_table_dict[[key]] <- c(drug_screen_table_dict[[key]], matched_wells[[well]]$reading)
        }
      }
    }
  }
  
  # Calculate the standard deviation using drug_screen_std
  drug_screen_table_std <- drug_screen_std(drug_screen_table_dict, row_index, col_index)
  
  # Create a data frame for the results using drug_combination_output
  result_data <- drug_combination_output(drug_screen_table_dict, drug_name2, drug_name1)
  
  return(list(result_data, drug_screen_table_std, drug_screen_table_dict))
}

filter_data_by_drug_combo <- function(
  plate_drug_printer_data, combo_wells, drug_combo
) {
  combo_wells_for_combo <- combo_wells[[toString(drug_combo)]]
  filtered_data <- subset(plate_drug_printer_data,
                          well %in% combo_wells_for_combo)
  return(filtered_data)
}