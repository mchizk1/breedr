#' A constructor function for breedr objects
#'
#' This function takes a data.frame containing breeding records and constructs a
#' breedr object for plotting pedigrees. Options are present for calculating
#' inbreeding statistics and coefficients of parentage.
#'
#' @param dataset A data.frame with exactly three columns: Individuals, female parents,
#' and male parents. The input column names are arbitrary and will be renamed in
#' the resulting breedr object.
#' @examples
#' habsburg_breedr <- new_breedr(habsburg)

new_breedr <- function(dataset){
  assertthat::assert_that(is.data.frame(dataset), msg = "input dataset must be a data.frame or something that may be coerced into one")
  assertthat::assert_that(ncol(dataset) == 3, msg = "Input data.frame must contain exactly three columns. See the habsburg dataset for an example.")
  colnames(dataset) <- c("Ind", "Par1", "Par2")
  dataset$Ind <- nameR(dataset$Ind)
  dataset$Par1 <- nameR(dataset$Par1)
  dataset$Par2 <- nameR(dataset$Par2)
  for(i in 1:nrow(dataset)){
    if (!(dataset$Par1[i] %in% dataset$Ind)){
      dataset <- rbind(dataset, c(dataset$Par1[i], NA, NA))
    }
    if (!(dataset$Par2[i] %in% dataset$Ind)){
      dataset <- rbind(dataset, c(dataset$Par2[i], NA, NA))
    }
  }
  dataset <- dataset[!is.na(dataset$Ind),]
  assertthat::assert_that(length(dataset$Ind) == length(unique(dataset$Ind)),
                          msg = paste0("Input data contains duplicate individuals: ",
                                       duplicated(dataset$Ind),
                                       collapse = " "))
  return(dataset)
}

