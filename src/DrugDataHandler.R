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

drug_combination_output <- function(
  drug_dict, row_drug, column_drug
) {
  sorted_drug_dict_keys <- sort_drug_dict_keys(drug_dict)

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
      paste0("Concentration ", row_drug, column_drug, sep = " "),
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
  colnames(final_data) <- c(paste0('Concentration ', row_drug, '_', column_drug), 
                            'Average_Reading', 'Median_Reading',
                            'STD_Reading', 'n', 'Raw Readings')
  
  return(final_data)
}

sort_drug_dict_keys <- function(
  drug_dict
) {
  drug_dict_keys <- names(drug_dict)
  tuples_keys <- lapply(drug_dict_keys, function(s) {
                              split_concentration <- as.numeric(
                                                      unlist(strsplit(s, "_")))
                              return(c(split_concentration[1], 
                                      split_concentration[2]))
                            })

  tuples_keys_matrix <- do.call(rbind, tuples_keys)

  order_indices <- order(tuples_keys_matrix[, 1])

  sorted_keys <- drug_dict_keys[order_indices]

  return(sorted_keys)
}

define_saved_data_file <- function(
  saved_path, analyze, plate, drug_names
  ) {
  file.path(saved_path, paste0(analyze,
                                gsub(" ", "_", plate), "_",
                                      drug_names[1], "_", drug_names[2],
                                      ".csv"  
                        )
  )
}
save_drug_data <- function(
  plate, matched_data, drug_data_table, output_dir
) {
  matched_wells <- matched_data[[1]]
  drug_conc_ranges <- matched_data[[2]]
  drug_screen_table_mean <- drug_data_table[[1]]
  drug_screen_table_std <- drug_data_table[[2]]
  drug_screen_table <- drug_data_table[[3]]


  drug_names <- names(drug_conc_ranges)
  saved_path <- file.path(output_dir, paste0(
                                          gsub(" ", "_", plate), "_", 
                                          drug_names[1], "_", 
                                          drug_names[2]))
  if (!file.exists(saved_path)) {
    dir.create(saved_path, showWarnings = FALSE)
  }
  saved_analyzed_file <- define_saved_data_file(saved_path, "Analyzed_",
                                                              plate, drug_names)
  write.csv(drug_screen_table_mean, saved_analyzed_file, row.names = FALSE, na = "0")
  saved_std_analyzed_file <- define_saved_data_file(saved_path,
                                                    "STD_Analyzed_",
                                                              plate, drug_names)
  write.csv(drug_screen_table_std, saved_std_analyzed_file,
                                                  row.names = FALSE, na = "0")
  saved_raw_analyzed_file <- define_saved_data_file(saved_path,
                                                    "Raw_Deconvoluted_",
                                                              plate, drug_names)
  drug_screen_table <- convert_match_matrix_2_frame(drug_screen_table)
  drug_screen_table <- format_match_data_frame(drug_screen_table)
  write.csv(drug_screen_table, saved_raw_analyzed_file,
                                                  row.names = FALSE, na = "0")
  
  return(drug_screen_table)
  #return(saved_path)
}

map_numbers_to_idx <- function(
  numbers_list
) {
  numbers_to_idx_dict <- list()
  for (idx in seq_along(numbers_list)) {
    numbers_to_idx_dict[as.character(numbers_list[idx])] <- idx
  }
  return(numbers_to_idx_dict)
}

create_data_table <- function(
  matched_data
) {
  matched_wells <- matched_data[[1]]
  drug_conc_ranges <- matched_data[[2]]
  #drug_names <- names(drug_conc_ranges)[!(names(drug_conc_ranges) %in% c("DMSO", "a+Tw", "a+B"))]
  drug_names <- names(drug_conc_ranges)
  if (length(drug_names) > 2) {
    if (!("DMSO" %in% drug_names) && !("a+Tw" %in% drug_names) && !("a+B" %in% drug_names)) {
      cat("Only two drug combinations can be processed, taking the first two drug names in:", drug_names[1:2], "\n")
    }
  }
  
  # Create tables for mean, n values, and standard deviation
  col_index <- as.numeric(drug_conc_ranges[[drug_names[1]]])
  row_index <- as.numeric(drug_conc_ranges[[drug_names[2]]])
  
  drug_screen_table <- matrix(0, nrow = length(row_index), ncol = length(col_index))
  drug_screen_table_n <- matrix(0, nrow = length(row_index), ncol = length(col_index))
  row_to_idx_dict <- map_numbers_to_idx(row_index)
  col_to_idx_dict <- map_numbers_to_idx(col_index)

  drug_screen_table_dict <- list()
  
  # Create enumerated drug_names dictionary for easy mapping
  drug_conc_ranges <- drug_conc_ranges[!(names(drug_conc_ranges) %in% c("DMSO", "a+Tw", "a+B"))]
  number_to_drug_name <- as.list(seq_along(drug_conc_ranges))
  drug_name_to_number <- setNames(number_to_drug_name, names(drug_conc_ranges))
  
  # Iterate over each well and drug_info
  for (well in names(matched_wells)) {
    drug_dictionary <- matched_wells[[well]]$drug_info
    ind <- 0.0
    col <- 0.0
    
    # Initialize values for situations when only one drug is added
    for (drug in names(drug_dictionary)) {
      if (!(drug %in% c("DMSO", "a+Tw", "a+B"))) {
        drug_location <- drug_name_to_number[[drug]]
        if (drug_location == 2) {
          ind <- drug_dictionary[[drug]]
        } else {
          col <- drug_dictionary[[drug]]
        }
      }
    }
    
    # Add value into proper cell in data tables
    matrix_ind <- row_to_idx_dict[[as.character(ind)]]
    matrix_col <- col_to_idx_dict[[as.character(col)]]

    drug_screen_table[matrix_ind, matrix_col] <- drug_screen_table[matrix_ind, matrix_col] + matched_wells[[well]]$reading
    drug_screen_table_n[matrix_ind, matrix_col] <- drug_screen_table_n[matrix_ind, matrix_col] + 1
    
    
    # Create easily split key for values for going back to the data table
    key <- paste0(as.character(ind), "_", as.character(col))
    drug_screen_table_dict[[key]] <- c(drug_screen_table_dict[[key]], matched_wells[[well]]$reading)
  }
  
  # Taking means and rounding to two decimal places
  drug_screen_table_mean <- round(drug_screen_table / drug_screen_table_n, digits = 2)
  drug_screen_table_std <- drug_screen_std(drug_screen_table_dict,
                                            row_index, col_index)
  drug_screen_table <- drug_combination_output(drug_screen_table_dict,
                                          drug_names[[2]], drug_names[[1]])

  return(list(drug_screen_table_mean, drug_screen_table_std, drug_screen_table))
}

convert_match_matrix_2_frame <- function(
  drug_data_matrix
) {
  rows <- rownames(drug_data_matrix)
  columns <- colnames(drug_data_matrix)
  num_columns <- length(columns)
  drug_data_frame <- data.frame(matrix(ncol = num_columns))
  colnames(drug_data_frame) <- columns
  for (row in rows) {
    record <- list()
    for (col in columns) {
        if (grepl("Concentration", col)) {
            record <- c(record, row)  
        } else {
            record <- c(record, drug_data_matrix[row, col])
        }
    }
    drug_data_frame <- rbind(drug_data_frame, record)
  }
  drug_data_frame <- na.omit(drug_data_frame)
  return(drug_data_frame)
}

format_match_data_frame <- function(
  drug_data_frame
) {
  max_n_reading = max(drug_data_frame$n)
  columns <- colnames(drug_data_frame)
  num_columns <- length(columns) - 1 + max_n_reading
  format_drug_data_frame <- data.frame(matrix(ncol = num_columns))
  format_columns_names <- columns[-length(columns)]
  for (idx in 1:max_n_reading) {
    col_name <- paste0("Raw Readings_", as.character(idx))
    format_columns_names <- c(format_columns_names, col_name)
  }
  colnames(format_drug_data_frame) <- format_columns_names
  for (row_idx in 1:nrow(drug_data_frame)){
    row_data <- drug_data_frame[row_idx, ]
    n_raw_readings <- row_data$n
    raw_readings <- row_data$`Raw Readings`
    raw_readings <- unlist(strsplit(raw_readings, ","))
    raw_readings <- as.numeric(raw_readings)
    record <- list()
    for (col in columns) {
        if (!grepl("Raw Readings", col)) {
          record <- c(record, row_data[[col]])
        }
    }
    for (raw_reading_idx in 1:max_n_reading) {
      if (raw_reading_idx < n_raw_readings) {
        record <- c(record, raw_readings[raw_reading_idx])
      } else {
        record <- c(record, 0)
      }
    }
    format_drug_data_frame <- rbind(format_drug_data_frame, record)
  }
  format_drug_data_frame <- format_drug_data_frame[-1, ]
  format_drug_data_frame <- sort_raw_deconvoluted_by_drug_two(
                                                        format_drug_data_frame)
  format_drug_data_frame <- sort_raw_deconvoluted_by_drug_one(
                                                          format_drug_data_frame)
  return(format_drug_data_frame)
}

