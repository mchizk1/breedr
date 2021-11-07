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
  }
  assertthat::assert_that(is.vector(genotypes), is.character(genotypes),
                          msg= "genotype must be a character vector or string")
  total <- 0
  for(i in 1:length(genotypes)){
    total = total + i
  }
  pb <- progress::progress_bar$new(format =
          "(:spin) [:bar] :percent [:current out of :total COP's calculated]",
                                   total = total,
                                   complete = "=",
                                   incomplete = "-",
                                   current = ">")
  COP <- matrix(nrow = length(genotypes), ncol = length(genotypes))
  genotypes <- nameR(genotypes, str_ops, na_val)
  rownames(COP) <- genotypes
  colnames(COP) <- genotypes
  diagonal <- 1
  for(i in genotypes){
    assertthat::assert_that(i %in% breedr$Ind,
                            msg = paste0(i, " was not found in the dataset provided."))
    idxi <- ped_idx(i, breedr = breedr)
    count <- 1
    COP_vec <- rep(0, length(genotypes))
    for(j in genotypes[1:diagonal]){
      pb$tick()
      idxj <- ped_idx(j, breedr = breedr)
      if(i %in% idxj$pedvec){
        COP_element <- sum(0.5^dplyr::filter(idxj, pedvec == i)$gidx)
      } else if (j %in% idxi$pedvec){
        COP_element <- sum(0.5^dplyr::filter(idxi, pedvec == j)$gidx)
      } else {
        COP_element <- 0
      }
      COP_vec[count] <- COP_element
      count = count+1
    }
    diagonal = diagonal+1
    COP[,i] <- COP_vec
  }
  COP <- COP + t(COP)
  for(i in 1:length(genotypes)){
    COP[i,i] = COP[i,i]/2
  }
  return(COP)
}


