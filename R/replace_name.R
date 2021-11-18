#' Name Conversions
#'
#' \code{replace_name()} replaces values by looking up name pairs from a key. This is
#' similar to Excel's vlookup function and may be helpful when switching from
#' selection codes to cultivar names.
#'
#' @param data_vec A character vector of ID codes to be replaced
#' @param key_old A character vector containing ID codes matching those of \code{data_vec}.
#' This vector may either be longer or shorter than \code{data_vec}.
#' @param key_new A character vector containing new ID codes corresponding to
#' replacement values for \code{key_old}. This vector must be of the same length
#' and order as \code{key_old}.
#' @return A character vector of length equal to \code{data_vec} containing newly
#' replaced ID codes
#' @examples
#' originals <- c("A-544", "A-730", "A-876", "A-1790")
#' key_old <- c("A-544", "A-730", "A-876")
#' key_new <- c("Cheyenne", "Shawnee", "Choctaw")
#' replace_name(originals, key_old, key_new)
#' @export

replace_name <- function(data_vec, key_old, key_new){
  assertthat::assert_that(is.character(data_vec), is.vector(data_vec),
              is.character(key_old), is.vector(key_old),
              is.character(key_new), is.vector(key_new),
              msg = "data_vec, key_old, and key_new must all be character vectors")
  assertthat::assert_that(length(key_old) == length(key_new), msg = "key_old and key_new must
              be of equal length")
  new_names <- c()
  for(i in 1:length(data_vec)){
    new_name <- data_vec[i]
    if(new_name %in% key_old){
      new_name = key_new[match(new_name, key_old)]
    }
    new_names <- c(new_names, new_name)
  }
  return(new_names)
}
