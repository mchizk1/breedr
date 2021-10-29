
#===============================================================================
# This function traces back a pedigree to the earliest available records and
# indexes the pedigree based on generation and position in generation
#===============================================================================

ped_idx <- function(parent, breedr){
  pedvec <- c(parent)
  counter <- 1
  gensz <-1
  nullcount <-0
  pidx <- c(1)
  gidx <- c(0)
  j <- 1
  repeat{
    k <- 0
    for(i in 1:gensz){
      g <- pedvec[counter]
      parents <- dplyr::filter(breedr, Ind == g)
      colnames(parents) <- NULL
      pedvec <- c(pedvec, unlist(parents[1,2:3]))
      k=k+2
      counter=counter+1
      if(is.na(parents[1,2]) | is.null(parents[1,2])){
        nullcount=nullcount+1
      }
      if(is.na(parents[1,3]) | is.null(parents[1,3])){
        nullcount=nullcount+1
      }
      pidx=c(pidx, k-1, k)
      gidx=c(gidx, j, j)
    }
    if(nullcount == 2*gensz){
      break
    }
    gensz=gensz*2
    nullcount=0
    j=j+1
  }
  #pedvec=pedvec[-1]
  ped_idx <- data.frame(pedvec=pedvec, pidx=pidx, gidx=gidx)
  ped_idx <- na.omit(ped_idx)
  return(ped_idx)
}

#===============================================================================
# This function takes two vectors of ancestors and returns unique matches
#===============================================================================

ca_detect <- function(P1_anc, P2_anc){
  ca <- c(NULL)
  for(i in P1_anc){
    if(i %in% P2_anc){
      ca = c(ca, i)
    }
  }
  return(ca)
}

#===============================================================================
# This function takes the results of two ped_idx() calls and identifies common
# ancestors between two individuals.  Then it uses index information to trace
# paths from one parent to the other, returning a vector of path lengths
# ("n" values).
#===============================================================================

ca_traceback <- function(P1, P2, breedr){
  P1_anc <- unique(P1$pedvec)
  P2_anc <- unique(P2$pedvec)
  # get common ancestors
  ca <- ca_detect(P1_anc, P2_anc)
  n_F <- c()
  if(!is.null(ca)){
    # Shrink vectors indices based on shared ancestry
    for(i in ca){
      P1 <- distinct_ca(P1, i, breedr)
      P2 <- distinct_ca(P2, i, breedr)
    }
    # Rerun common ancestors after removals
    P1_anc <- unique(P1$pedvec)
    P2_anc <- unique(P2$pedvec)
    ca <- ca_detect(P1_anc, P2_anc)
    traces <- c()
    F_vec <- c()
    for(i in ca){
      P1_scan <- dplyr::filter(P1, pedvec == i)
      P2_scan <- dplyr::filter(P2, pedvec == i)
      F_stat <- dplyr::filter(breedr, Ind == i)[1,4]
      # For each occurrence on the P1 tree...
      for(j in 1:nrow(P1_scan)){
        # For each occurrence on the P2 tree...
        for(k in 1:nrow(P2_scan)){
          trace_n <- P1_scan$gidx[j] + P2_scan$gidx[k] + 1
          traces = c(traces, trace_n)
          F_vec = c(F_vec, F_stat)
        }
      }
    }
    n_F <- data.frame(traces=traces, F_stat=F_stat)
  }
  return(n_F)
}

#===============================================================================
# When a common ancestor is identified, all preceding ancestors in ped_idx
# will need to be filtered out to avoid inflation of F with redundant
# common ancestry information.
#===============================================================================

distinct_ca <- function(Par_idx, ca, breedr){
  if(ca %in% Par_idx$pedvec){
    ca_group <- dplyr::filter(Par_idx, pedvec == ca)
    # For each occurrence of common ancestor i...
    for(i in 1:nrow(ca_group)){
      # If the common ancestor is not a dead end...
      if (ca_group$gidx[i] < max(Par_idx$gidx)){
        # For each generation behind common ancestor...
        for(j in 1:(max(Par_idx$gidx) - ca_group$gidx[i])){
          gen_up <- ca_group$gidx[i]+j
          parnum <- ca_group$pidx[i]*2^j
          Par_idx <- dplyr::filter(Par_idx, !((gidx == gen_up) &
                                                (pidx %in% c((parnum-2^j+1):(parnum)))))
        }
      }
    }
  }
  return(Par_idx)
}

#===============================================================================
# This function calls ped_idx() and ca_traceback() to assemble necessary data
# and calculates Wright's F stat for an individual in the breedr object
#===============================================================================

#' Wright's F statistic
#'
#' This function calculates Wright's inbreeding F statistic for any set of
#' parentage data. It also calls \code{NameR()} to standardize name input for flexibility.
#'
#' @param breedr A data.frame containing three columns: individual, female parent,
#' and male parent
#' @return A numeric scalar or vector (Wright's F) for the requested genotype(s)

wrights_f <- function(breedr, iteration){
  genotype = breedr$Ind
  genotype <- nameR(genotype)
  assertthat::assert_that(is.vector(genotype), is.character(genotype),
                          msg= "genotype must be a character vector or string")
  F_unadjusted <- c()
  pb <- progress::progress_bar$new(format =
              "(:spin) [:bar] :percent [Step :it out of 2: :current out of :total F's calculated]",
                         total = length(genotype),
                         complete = "=",
                         incomplete = "-",
                         current = ">")
  for(i in genotype){
    assertthat::assert_that(i %in% breedr$Ind,
                            msg = paste0(i, " was not found in the breedr object
                                         provided"))
    pb$tick(tokens = list("it" = iteration))
    parents <- dplyr::filter(breedr, Ind == i)
    # This conditional handles selfing (same parent in both columns)
    if(!is.na(parents[1,2]) & !is.na(parents[1,3]) &
       parents[1,2] == parents[1,3]){
      n_F <- data.frame(traces=1, F_stat=dplyr::filter(breedr,
                                                    Ind == parents[1,2])[1,4])
    } else {
      P1_ped <- ped_idx(parents$Par1, breedr)
      P2_ped <- ped_idx(parents$Par2, breedr)
      n_F <- ca_traceback(P1_ped, P2_ped, breedr)
    }
    n_F$F_stat <- as.numeric(n_F$F_stat)
    F_contribution <- 0.5^(n_F$traces)*(1+n_F$F_stat)
    F_i <- sum(F_contribution)
    F_unadjusted <- c(F_unadjusted, F_i)
  }
  return(F_unadjusted)
}

#===============================================================================
# This function calls wrights_f to calculate f for all individuals in a dataset
# in a two step procedure
#===============================================================================

#' Wright's F statistic
#'
#' Using pedigree data, this function calculates Wright's
#' inbreeding F statistic for a set of parentage data. It also calls
#' \code{NameR()} to standardize name input for flexibility.
#'
#' @param breedr A data.frame containing three columns: individual, female parent,
#' and male parent
#' @return A data.frame containing inbreeding statistics of all individuals in
#' the input data.frame
#' @examples
#' breedr_F(habsburg)
#' @export

breedr_F <- function(breedr){
  breedr <- new_breedr(breedr)
  breedr$F_stat <- 0
  breedr$F_stat <- wrights_f(breedr, 1)
  breedr$F_stat <- wrights_f(breedr, 2)
  return(breedr)
}
