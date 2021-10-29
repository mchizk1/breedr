#' Standardizing Names.
#'
#' \code{nameR} takes potentially inconsistent or erroneous names of
#' cultivars and breeding selections, and replaces them with standardized
#' alternatives. To minimize mismatched data or dropped entries, one could run
#' all occurrences of names through \code{nameR} prior to analysis.
#'
#' @param names Either a single string or a vector of character strings containing the
#' names of cultivars and/or breeding selections.
#' @return Either a single string or a vector of character strings containing the
#' standardized name replacements. Length will be equivalent to the length of
#' \code{names}.
#'
#' @export

nameR <- function(names, strip_ws = T, cap_sns = F) {
  assertthat::assert_that(is.character(names), msg = "names must be a
                          character string or vector of character strings")
  if(strip_ws == T){
    # Strip leading/trailing white space
    names <- stringr::str_trim(names, side = "both") %>%
      # Remove repeated internal white space
      stringr::str_squish()
  }
  if(cap_sns == F){
    # Make all uppercase
    names <- toupper(names)
  }
  # Replace all spaces with underscores
  names <- stringr::str_replace_all(names, "[:blank:]+", "_")
  return(names)
}
