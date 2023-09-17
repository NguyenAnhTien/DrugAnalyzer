# @date: 2023-09-17

split_string <- function(
  string
) {
  parts <- unlist(strsplit(string, " "))  # Split the string by spaces
  if (length(parts) > 0) {
    return(parts[1])  # Return the first element
  } else {
    return(string)  # If no splitting occurred, return the original value
  }
}
