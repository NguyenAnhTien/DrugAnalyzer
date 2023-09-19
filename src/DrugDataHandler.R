# @date: 2023-09-18

cal_drug_combinations <- function(
  drug_printer_data
) {
  # @desc:
  # - Only handle two drugs
  # - However, there are some duplication
  # - The lines of code below the line 16 are used to handle the duplication
  plate_drugs <- unique(drug_printer_data$drug)
  plate_drugs_cleaned <- plate_drugs[
                                     !(plate_drugs %in% c(
                                                          "DMSO",
                                                          "a+Tw", "a+B"))]
  solvent <- plate_drugs[plate_drugs %in% c("DMSO", "a+Tw", "a+B")]

  if (length(plate_drugs_cleaned) == 2) {
    return(list(c(plate_drugs_cleaned, solvent)))
  }

  combination_wells <- unique(drug_printer_data$well[
                                                     duplicated(
                                                      drug_printer_data$well)])

  drug_printer_data_ordered <- drug_printer_data[order(
                                              drug_printer_data$drug), ]

  drug_combinations <- list()
  for (well in combination_wells) {
    well_data <- subset(drug_printer_data_ordered, well == well)
    well_drugs <- unique(well_data$drug)
    if (!identical(well_drugs, drug_combinations)) {
      drug_combinations <- append(drug_combinations, list(well_drugs))
    }
  }

  drug_combinations <- drug_combinations[sapply(drug_combinations, length) > 2]

  return(drug_combinations)
}

cal_drug_combo_wells <- function(drug_printer_data, drug_combinations) {
  drug_combo_wells <- list()
  for (drug_combo in drug_combinations) {
    all_wells_with_drug <- unique(
                                  drug_printer_data$well[
                                      drug_printer_data$drug %in% drug_combo])
    drug_combo_wells[[toString(drug_combo)]] <- list()
    
    for (well in all_wells_with_drug) {
      well_drugs <- drug_printer_data$drug[drug_printer_data$well == well]
      well_drugs <- as.character(well_drugs)
      correct_well <- TRUE
      
      for (drug in well_drugs) {
        if (!(drug %in% drug_combo)) {
          correct_well <- FALSE
          break
        }
      }
      
      if (correct_well) {
        drug_combo_wells[[
                          toString(drug_combo)]] <- append(
                                                            drug_combo_wells[[
                                                            toString(
                                                            drug_combo)]], well)
      }
    }
  }
  
  return(drug_combo_wells)
}

drug_screen_std <- function(drug_screen_dict, row_index, col_index) {
  # Create an empty data frame with specified row and column indices
  drug_screen_standard_dev <- matrix(0, 
                                     nrow = length(row_index),
                                     ncol = length(col_index))
  rownames(drug_screen_standard_dev) <- row_index
  colnames(drug_screen_standard_dev) <- col_index
  
  for (drug_combination in names(drug_screen_dict)) {
    standard_dev <- sd(unlist(drug_screen_dict[[drug_combination]]), na.rm = TRUE)
    split_concentration <- as.numeric(unlist(strsplit(drug_combination, "_")))
    ind <- split_concentration[1]
    col <- split_concentration[2]
    
    if (!is.na(standard_dev)) {
      drug_screen_standard_dev[ind, col] <- standard_dev
    } else {
      drug_screen_standard_dev[ind, col] <- 0
    }
  }
  
  return(as.data.frame(drug_screen_standard_dev))
}

drug_combination_output <- function(drug_dict, row_drug, column_drug) {
  # Sort the drug dictionary keys
  sorted_drug_dict_keys <- sort(names(drug_dict), 
                                 index.return = TRUE, 
                                 fun = function(x) {
                                   split_concentration <- as.numeric(
                                                      unlist(strsplit(x, "_")))
                                   return(c(split_concentration[1], 
                                            split_concentration[2]))
                                 })$x

  # Initialize a list to store the results
  average_readings <- list()
  
  for (drug_pair in sorted_drug_dict_keys) {
    readings <- unlist(drug_dict[[drug_pair]])
    mean_reading <- mean(readings, na.rm = TRUE)
    median_reading <- median(readings, na.rm = TRUE)
    std_reading <- sd(readings, na.rm = TRUE)
    n <- length(readings)
    raw_readings <- paste(readings, collapse = ', ')
    
    # Create a list for this drug pair
    drug_pair_data <- list(
      paste("Concentration", row_drug, column_drug, sep = " "),
      mean_reading,
      median_reading,
      std_reading,
      n,
      raw_readings
    )
    
    # Append the drug pair data to the results list
    average_readings[[drug_pair]] <- drug_pair_data
  }

  # Convert the list to a data frame
  final_data <- do.call(rbind, average_readings)
  colnames(final_data) <- c(f'Concentration {row_drug}_{column_drug}', 
                            'Average_Reading', 'Median_Reading', 
                            'STD_Reading', 'n', 'Raw Readings')
  
  return(final_data)
}
