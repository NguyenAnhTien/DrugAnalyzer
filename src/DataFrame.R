# @date: 2023-09-17
source("Utils.R")

rename_columns <- function(
  df, column_names_mapping
) {
  rename_df <- setNames(df, column_names_mapping)
  return(rename_df)
}

sort_raw_deconvoluted_by_drug_two <- function(df) {
  column <- names(df)[1]
  data <- list()
  
  for (idx in 1:nrow(df)) {
    row <- df[idx, , drop = FALSE]
    value <- as.numeric(unlist(strsplit(row[[column]], "_"))[2])
    data[[idx]] <- list(record = row, value = value)
  }
  
  values <- sapply(data, function(item) item$value)
  sorted_indices <- order(values)
  data <- data[sorted_indices]
  
  df_ <- do.call(rbind, lapply(data, function(x) x$record))
  
  return(df_)
}

sort_raw_deconvoluted_by_drug_one <- function(df) {
  column <- names(df)[1]
  columns <- names(df)
  data <- list()
  
  for (idx in 1:nrow(df)) {
    row <- df[idx, , drop = FALSE]
    values <- unlist(strsplit(row[[column]], "_"))
    drug_one <- values[1]
    drug_two <- values[2]
    
    if (!is.na(drug_one) && !is.na(drug_two)) {
        if (drug_one %in% names(data)) {
            data[[drug_one]] <- append(data[[drug_one]], list(c(record = row,
                                                value = as.numeric(drug_two))))
        } else {
            data[[drug_one]] <- list(c(record = row,
                                                value = as.numeric(drug_two)))
        }
    }
  }
  
  records <- list()
  for (key in names(data)) {
    group_items <- data[[key]]
    sorted_group_items <- sort_dict_values(group_items)
    for (item in sorted_group_items) {
        record <- item[-which(names(item) == "value")]
        records <- append(records, list(record))
    }
  }
  
  df_ <- do.call(rbind, records)
  colnames(df_) <- columns
  return(df_)
}
