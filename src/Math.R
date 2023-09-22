# @date: 2023-09-17
# @desc:
#   - This module contains functions for calculating approximate numbers

round_sig <- function(
  x, sig = 1, small_value = 1.0e-9
) {
  tryCatch(
    {
      # Calculate the number of significant figures
      sig_fig <- sig - floor(log10(max(abs(as.numeric(x)),
                                       abs(small_value)))) - 1

      # Round to the specified number of significant figures
      return(round(as.numeric(x), sig_fig))
    },
    error = function(e) {
      return(x)  # Return the original value if an error occurs
    }
  )
}