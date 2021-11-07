#' Standardizing Names.
#'
#' \code{nameR} takes potentially inconsistent or erroneous names of
#' cultivars and breeding selections, and replaces them with standardized
#' alternatives. To minimize mismatched data or dropped entries, one could run
#' all occurrences of names through \code{nameR} prior to analysis.
#'
#' @param names Either a single string or a vector of character strings containing the
#' names of cultivars and/or breeding selections.
#' @param str_ops A character vector containing optional methods for string manipulation.
#' Possible options include: \cr
#' 'strip_ws' (default) - Strip leading and trailing white space. \cr
#' 'upper' (default) - Make all uppercase to eliminate case differences.\cr
#' 'lower' - Make all lowercase to eliminate case differences.  Cannot be combined with 'upper'. \cr
#' 'shrink_ws' (default) - remove duplicated white space from center of strings. \cr
#' NULL is acceptable to skip string manipulation
#' @param na_val An optional string indicating a value to be treated as missing.
#' @return Either a single string or a vector of character strings containing the
#' standardized name replacements. Length will be equivalent to the length of
#' \code{names}.
#'
#' @export

nameR <- function(names, str_ops = c("strip_ws", "upper", "shrink_ws"), na_val = NULL) {
  avail_ops <- c("upper", "lower", "strip_ws", "shrink_ws")
  assertthat::assert_that(is.character(names), msg = "names must be a
                          character string or vector of character strings")
  if(!is.null(str_ops)){
    assertthat::assert_that(is.vector(str_ops) & is.character(str_ops),
                            msg = "str_ops must be a character vector")
    str_logical <- str_ops %in% avail_ops
    assertthat::assert_that(!(F %in% str_logical),
                            msg = "available options for str_ops include 'upper', 'lower', 'strip_ws', and 'shrink_ws'")
    assertthat::assert_that(!("upper" %in% str_ops & "lower" %in% str_ops),
                            msg = "str_ops cannot include both 'upper' and 'lower'")
    if(!is.null(na_val)){
      assertthat::assert_that(assertthat::is.string(na_val), msg = "na_val must be a string")
      names[names == na_val] <- NA
    }
    if("strip_ws" %in% str_ops){
      # Strip leading/trailing white space
      names <- stringr::str_trim(names, side = "both")
    }
    if("shrink_ws" %in% str_ops){
      # Remove repeated internal white space
      names <- stringr::str_squish(names)
    }
    if("upper" %in% str_ops){
      # Make all uppercase
      names <- toupper(names)
    }
    if("lower" %in% str_ops){
      # Make all uppercase
      names <- tolower(names)
    }
  }
  # Replace all spaces with underscores
  names <- stringr::str_replace_all(names, "[:blank:]", "_")
  return(names)
}
