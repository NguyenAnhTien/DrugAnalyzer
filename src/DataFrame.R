# @date: 2023-09-17

rename_columns <- function(
  df, column_names_mapping
) {
  rename_df <- setNames(df, column_names_mapping)
  return(rename_df)
}