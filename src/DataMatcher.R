
match_data <- function(
  drug_printer_data, plate_reader_data
) {
  matched_wells <- list()

  unique_wells <- unique(drug_printer_data$well)

  for (well in unique_wells) {
    well_data <- subset(drug_printer_data, well == well)
    matched_well <- list(
      drug_info = as.list(
                          setNames(well_data$concentration, well_data$drug)),
      reading = as.integer(plate_reader_data$result[
                                                plate_reader_data$well == well])
    )

    matched_wells[[well]] <- matched_well
  }
  if (length(plate_reader_data$well) == 96) {
    letters <- c("B", "C", "D", "E", "F", "G")
    numbers <- c("02", "03", "04", "05", "06", "07", "08", "09", "10", "11")
    all_wells <- outer(letters, numbers, FUN = paste0)

    control_well <- NULL

    difference <- setdiff(all_wells, unique(drug_printer_data$well))

    if (length(difference) > 0) {
      control_well <- difference[1]
      matched_wells[[control_well]] <- list(
        drug_info = as.list(setNames(
                                     rep(0.0,
                                         length(
                                           unique(drug_printer_data$drug)
                                         )),
                                     unique(drug_printer_data$drug))),
        reading = as.integer(
                             plate_reader_data$result[
                                        plate_reader_data$well == control_well])
      )
    }
  }
  drug_concentration_ranges <- list()

  unique_drugs <- unique(drug_printer_data$drug)

  for (drug_name in unique_drugs) {
    drug_specific_concentrations <- unique(
                                           drug_printer_data$concentration[
                                        drug_printer_data$drug == drug_name])

    sorted_concentrations <- sort(c(drug_specific_concentrations, 0))

    drug_concentration_ranges[[drug_name]] <- sorted_concentrations
  }
  return(list(matched_wells, drug_concentration_ranges))
}