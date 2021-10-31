#' Coefficients of Parentage for UA blackberry program
#'
#' Using pedigree data, this function computes a matrix containing
#' coefficients of parentage (COP) for any group of cultivars or breeding selections.
#' It also calls \code{NameR()} to standardize name input for flexibility.
#'
#' @param breedr A data.frame containing three columns: individual, female parent,
#' and male parent
#' @return A numeric square matrix containing COP's for the requested
#' @examples
#' breedr_COP(habsburg)
#' @export

breedr_COP <- function(breedr){
  breedr <- new_breedr(breedr)
  genotypes = breedr$Ind
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
  genotypes <- nameR(genotypes)
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


