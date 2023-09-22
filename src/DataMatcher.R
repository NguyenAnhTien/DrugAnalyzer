# @date: 2023-09-22
# @desc: This module contains functions for comparing drug data

match_data <- function(
  drug_printer_data, plate_reader_data
) {
  matched_wells <- list()

  unique_wells <- unique(drug_printer_data$well)

  for (well in unique_wells) {
    well_data <- drug_printer_data[drug_printer_data$well == well, ]
    matched_well <- list(
      drug_info = as.list(
                          setNames(well_data$concentration, well_data$drug)),
      reading = as.integer(plate_reader_data$result[
                                                plate_reader_data$well == well])
    )

    matched_wells[[well]] <- matched_well
  }
  
  drug_concentration_ranges <- list()

  unique_drugs <- unique(drug_printer_data$drug)

  for (drug_name in unique_drugs) {
    if (!(drug_name %in% c("DMSO", "a+Tw", "a+B"))){
      drug_specific_concentrations <- unique(
                                           drug_printer_data[
                                            drug_printer_data$drug == drug_name,
                                        ]$concentration)

      sorted_concentrations <- sort(c(drug_specific_concentrations, 0))

      drug_concentration_ranges[[drug_name]] <- sorted_concentrations
    }
  }
  return(list(matched_wells, drug_concentration_ranges))
}