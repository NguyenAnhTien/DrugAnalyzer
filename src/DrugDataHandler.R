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
  colnames(final_data) <- c(paste('Concentration', row_drug, '_', column_drug), 
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
                                    substr(
                                            unlist(
                                              strsplit(drug_names[1], " ")
                                            ), 1, 1
                                    ), "_", 
                                    substr(
                                            unlist(
                                              strsplit(drug_names[2], " ")
                                            ), 1, 1
                                    ),
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
  print(drug_names)
  saved_path <- file.path(output_dir, paste0(
                                          gsub(" ", "_", plate), "_", 
                                          substr(
                                            unlist(
                                              strsplit(drug_names[1], " ")
                                            ), 1, 1
                                          ), "_", 
                                          substr(
                                            unlist(
                                              strsplit(drug_names[2], " ")
                                            ), 1, 1)))
  if (!file.exists(saved_path)) {
    dir.create(saved_path, showWarnings = FALSE)
  }
  saved_analyzed_file <- define_saved_data_file(saved_path, "Analyzed_",
                                                              plate, drug_names)
  write.csv(drug_screen_table_mean, saved_analyzed_file, row.names = TRUE, na = "0")
  saved_std_analyzed_file <- define_saved_data_file(saved_path,
                                                    "STD_Analyzed_",
                                                              plate, drug_names)
  write.csv(drug_screen_table_std, saved_std_analyzed_file,
                                                  row.names = TRUE, na = "0")
  saved_raw_analyzed_file <- define_saved_data_file(saved_path,
                                                    "Raw_Deconvoluted_",
                                                              plate, drug_names)
  write.csv(drug_screen_table, saved_raw_analyzed_file,
                                                  row.names = TRUE, na = "0")
}

drug_screen_viability <- function(
  drug_data_table, output_path, drug_printer_file_path, gr = FALSE
) {
  drug_screen_table <- drug_data_table[[3]]
  control_reading <- drug_screen_table[1, 1]
  
  if (gr) {
    viability <- drug_screen_table
  } else {
    viability <- (drug_screen_table / control_reading) * 100
  }
  
  # Read concentration unit
  file_data <- read.xlsx(drug_printer_file_path, sheetName = "Tabular detail")
  conc_unit <- unique(na.omit(file_data$`Conc. units`))
  
  if (length(conc_unit) > 1) {
    stop("Only one drug concentration unit can be used for SynergyFinder.")
  }
  
  write.csv(viability, file = output_path, row.names = FALSE, na = 0)
  
  if (conc_unit[1] != "nM") {
    # Convert concentration to nM for proper axis labeling in SynergyFinder
    viability_nm <- viability
    
    # Scale concentration to nM
    colnames(viability_nm) <- colnames(viability_nm) * 1000
    rownames(viability_nm) <- rownames(viability_nm) * 1000
    
    # Write to CSV with nM unit
    output_path_nm <- paste0(tools::file_path_sans_ext(output_path),
                                                            "_nM_Converted.csv")
    write.csv(viability_nm, file = output_path_nm, row.names = FALSE, na = 0)
  }
  
  return(viability)
}

create_data_table <- function(
  matched_data
) {
  matched_wells <- matched_data[[1]]
  drug_conc_ranges <- matched_data[[2]]
   # Extract drug names excluding DMSO, a+Tw, and a+B
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
  rownames(drug_screen_table) <- row_index
  colnames(drug_screen_table) <- col_index
  
  drug_screen_table_n <- matrix(0, nrow = length(row_index), ncol = length(col_index))
  rownames(drug_screen_table_n) <- row_index
  colnames(drug_screen_table_n) <- col_index
  
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
        
        if (drug_location) {
          ind <- drug_dictionary[[drug]]
        } else {
          col <- drug_dictionary[[drug]]
        }
      }
    }
    
    # Add value into proper cell in data tables
    tryCatch({
      drug_screen_table[ind, col] <- drug_screen_table[ind, col] + matched_wells[[well]]$reading
      drug_screen_table_n[ind, col] <- drug_screen_table_n[ind, col] + 1
    }, error = function(e) {
      cat("\nCheck to make sure single or multi-parameter is correct.\n")
      stop(e)
    })
    
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
  
  # drug_screen_table_mean <- as.data.frame(drug_screen_table_mean)
  # drug_screen_table_std <- as.data.frame(drug_screen_table_std)
  # drug_screen_table <- as.data.frame(drug_screen_table)

  return(list(drug_screen_table_mean, drug_screen_table_std, drug_screen_table))
}

format_drug_data <- function(
  drug_data_table
) {
  rows <- rownames(drug_data_table)
  columns <- colnames(drug_data_table)
  num_columns <- length(columns)
  drug_data_frame <- data.frame(matrix(ncol = num_columns))
  colnames(drug_data_frame) <- columns
  for (row in rows) {
    record <- list()
    for (col in columns) {
        if (grepl("Concentration", col)) {
            record <- c(record, row)  
        } else {
            record <- c(record, drug_data_table[row, col])
        }
    }
    drug_data_frame <- rbind(drug_data_frame, record)
  }
  return(drug_data_frame)
}