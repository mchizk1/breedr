#' Coefficients of Parentage for UA blackberry program
#'
#' Using pedigree data, this function computes a matrix containing
#' coefficients of parentage (COP) for any group of cultivars or breeding selections.
#' It also calls \code{NameR()} to standardize name input for flexibility.
#'
#' @param breedr A data.frame containing three columns: individual, female parent,
#' and male parent
#' @param genotypes A vector of character strings indicating a subset of individuals
#' to calculate COPs for. By default, all individuals in the breedr dataset are included.
#' @param str_ops A character vector containing optional methods for string manipulation.
#' Possible options include: \cr
#' 'strip_ws' (default) - Strip leading and trailing white space. \cr
#' 'upper' (default) - Make all uppercase to eliminate case differences.\cr
#' 'lower' - Make all lowercase to eliminate case differences.  Cannot be combined with 'upper'. \cr
#' 'shrink_ws' (default) - remove duplicated white space from center of strings. \cr
#' NULL is acceptable to skip string manipulation
#' @param na_val An optional string indicating a value to be treated as missing.
#' @return A numeric square matrix containing COP's for the requested
#' @examples
#' breedr_COP(habsburg)
#' @export

breedr_COP <- function(breedr, genotypes = NULL,
                       str_ops = c("strip_ws", "upper", "shrink_ws"), na_val = NULL){
  breedr <- new_breedr(breedr, str_ops, na_val)
  if(is.null(genotypes)){
    genotypes = breedr$Ind
  } else {
    genotypes <- nameR(genotypes, str_ops = str_ops, na_val = na_val)
  }
  assertthat::assert_that(is.vector(genotypes), is.character(genotypes),
                          msg= "genotype must be a character vector or string")
  idx <- purrr::map(genotypes, ped_idx, breedr = breedr)
  names(idx) <- genotypes
  cross_list <- list(genotypes, genotypes)
  names(cross_list) <- c("a", "b")
  pairs <- purrr::cross_df(cross_list)
  COP <- purrr::map2(.x = pairs$a, .y = pairs$b, .f = COP_compute, idx = idx) %>%
    unlist() %>%
    matrix(nrow = length(genotypes))
  rownames(COP) <- genotypes
  colnames(COP) <- genotypes
  return(COP)
}

COP_compute <- function(a, b, idx){
  if(b %in% idx[[a]]$pedvec){
    COP_element <- sum(0.5^dplyr::filter(idx[[a]], pedvec == b)$gidx)
  } else if (a %in% idx[[b]]$pedvec){
    COP_element <- sum(0.5^dplyr::filter(idx[[b]], pedvec == a)$gidx)
  } else {
    COP_element <- 0
  }
  return(COP_element)
}
