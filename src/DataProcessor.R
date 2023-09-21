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
      return(c(match_data=matched_data, drug_data_table=drug_data_table))
      #save_drug_data(plate, matched_data, drug_data_table, output_dir)
    #   viability_path <- file.path(save_path, paste("Viability",
    #                                               gsub(" ", "_", plate), 
    #                               gsub(" ", "_", substr(drug_names[1], 1,
    #                                                 regexpr(" ", 
    #                                                 drug_names[1])-1)), ".csv"))

    #   viability <- drug_screen_viability(drug_data_table, viability_path,
    #                                                         drug_printer_file)
    #   synergy_finder_r_package(viability, drug_printer_file,
    #                           plate <- gsub(" ", "_", plate), saved_path)
    #   viability_tables <- c(viability_tables, list(viability))
    # synergyfinder(viability_tables, output_dir, drug_printer_file)
    }
  }
}

synergy_finder_r_package <- function(
  viability_table, drug_printer_file_path, plate, output_directory_path
) {
  # Getting concentration unit of liquid
  file_data <- read.xlsx(drug_printer_file_path, sheetName = "Tabular detail")
  conc_unit <- unique(na.omit(file_data$`Conc. units`))
  
  if (length(conc_unit) > 1) {
    stop("Only one drug concentration unit can be used for SynergyFinder.")
  }
  
  # Getting drug printer name for naming of synergy finder folder
  drug_printer_name <- tools::file_path_sans_ext(basename(drug_printer_file_path))
  
  drug1 <- rownames(viability_table)
  drug2 <- colnames(viability_table)
  
  synergy_finder_r_file <- file.path(output_directory_path, 
                                      paste0("SynergyFinder_",
                                      gsub(" ", "_", drug_printer_name),
                                      "_", drug1, "_", drug2, "_",
                                      plate, ".csv"))
  
  # Stack the table
  stacked_table <- as.data.frame.table(viability_table)
  colnames(stacked_table) <- c("ConcRow", "ConcCol", "Response")
  stacked_table$BlockId <- plate
  stacked_table$DrugRow <- drug1
  stacked_table$DrugCol <- drug2
  stacked_table$conc_r_unit <- conc_unit
  stacked_table$conc_c_unit <- conc_unit
  
  # Reorder columns
  stacked_table <- stacked_table[, c("BlockId", "DrugRow", "DrugCol",
                                    "ConcRow", "ConcCol", "Response",
                                    "conc_r_unit", "conc_c_unit")]
  
  # Write to CSV
  write.csv(stacked_table, file = synergy_finder_r_file,
            row.names = FALSE, sep = ",")
  
  return(stacked_table)
}

synergyfinder <- function(viability_tables, output_directory_path, drug_printer_file_path) {
  # Getting drug printer name for naming of synergy finder folder
  drug_printer_name <- tools::file_path_sans_ext(basename(drug_printer_file_path))

  # Getting concentration unit of liquid
  file_data <- read.xlsx(drug_printer_file_path, sheetName = "Tabular detail")
  conc_unit <- unique(na.omit(file_data$`Conc. units`))
  
  if (length(conc_unit) > 1) {
    stop("Only one drug concentration unit can be used for SynergyFinder.")
  }

  synergy_finder_matrix_file <- file.path(output_directory_path, 
                                           paste0("SynergyFinder_", gsub(" ", "_", drug_printer_name), ".csv"))

  # Create a list to hold the CSV data
  csv_data <- list()
  
  for (v_tab in viability_tables) {
    drug1 <- rownames(v_tab)
    drug2 <- colnames(v_tab)
    
    # Create a CSV data frame for this table
    csv_df <- data.frame(matrix(NA, nrow = nrow(v_tab) + 2, ncol = ncol(v_tab)))
    rownames(csv_df) <- NULL
    
    # Fill in the CSV data frame
    csv_df[1, ] <- c("Drug1:", drug1)
    csv_df[2, ] <- c("Drug2:", drug2)
    csv_df[3, ] <- c("ConcUnit:", conc_unit)
    csv_df[4:(nrow(v_tab) + 3), ] <- v_tab
    
    csv_data <- append(csv_data, list(csv_df))
  }
  
  # Write the list of CSV data frames to a single CSV file
  write.csv(do.call(rbind, csv_data), file = synergy_finder_matrix_file, row.names = FALSE)
}


filter_data_by_drug_combo <- function(
  plate_drug_printer_data, combo_wells, drug_combo
) {
  combo_wells_for_combo <- combo_wells[[toString(drug_combo)]]
  filtered_data <- subset(plate_drug_printer_data,
                          well %in% combo_wells_for_combo)
  return(filtered_data)
}