
gamete_ratio <- function(parent, ploidy){
  stringr::str_split(parent, "") %>%
    unlist() %>%
    combn(ploidy/2) %>%
    as.data.frame()
}

#' Expected segregation ratios.
#'
#' This function takes parental genotypes and projects expected segregation ratios of the progeny.
#'
#' @param parent1 A character string or positive integer indicating the genotype of the first parent.  If
#' an integer is provided, ploidy must be provided as well.  If a character format is used, ploidy will
#' be inferred. Only biallelic SNPs are currently supported.
#' @param parent2 A character string or positive integer indicating the genotype of the second parent.  Must be of
#' the same class as parent1. If not provided, assumed to be a self of parent1.
#' @param ploidy A positive integer indicating the ploidy level of the parents and progeny.  Only required if
#' parental genotypes are integers.
#' @return A data.frame providing expected segregation ratios of progeny.
#' @export
#' @examples
#' seg_ratio("Aa", "AA")
#' seg_ratio(1, 2, 4)
#' seg_ratio("AAAaaa")

seg_ratio <- function(parent1, parent2=parent1, ploidy=NULL){
  label <- parent1
  if((is.numeric(parent1) & is.numeric(parent2)) |
     (is.integer(parent1) & is.integer(parent2))){
    assertthat::assert_that(!is.null(ploidy),
                            msg = "Ploidy level must be specified if numeric data is used for parental genotypes")
    parent1 <- stringr::str_flatten(c(rep("a", parent1),
                             rep("b", ploidy-parent1)))
    parent2 <- stringr::str_flatten(c(rep("a", parent2),
                             rep("b", ploidy-parent2)))
    num_type <- T
  } else if(is.character(parent1) & is.character(parent2)){
    assertthat::assert_that(length(unlist(stringr::str_split(parent1, ""))) == length(unlist(stringr::str_split(parent2, ""))),
                            msg = "Parents must share the same ploidy level.")
    ploidy = length(unlist(stringr::str_split(parent1, "")))
    num_type <- F
  }
  alleles <- c(unlist(stringr::str_split(parent1, "")),
               unlist(stringr::str_split(parent2, ""))) %>%
    unique() %>%
    sort(decreasing = T)
  assertthat::assert_that(length(alleles) <= 2, msg = "Locus cannot have more than two alleles")
  if(length(alleles)==2){
    gametes1 <- gamete_ratio(parent1, ploidy)
    gametes2 <- gamete_ratio(parent2, ploidy)
    classes <- data.frame(Genotype = unlist(purrr::map(as.list(0:ploidy), ~ paste0(
      stringr::str_flatten(rep(alleles[1], ploidy-.)),
      stringr::str_flatten(rep(alleles[2], .))))
    ))
    geno_mat <- purrr::cross2(gametes1, gametes2) %>%
      purrr::map(~ unlist(.)  %>%
            sort(decreasing = T) %>%
            stringr::str_flatten()) %>%
      unlist() %>%
      table %>%
      as.data.frame() %>%
      dplyr::mutate(Freq = reduce_gcd(Freq)) %>%
      dplyr::rename(Genotype = '.') %>%
      dplyr::full_join(data.frame(Genotype = classes),
                by = "Genotype")
    if(num_type){
      geno_mat$Genotype <- stringr::str_count(geno_mat$Genotype, "a")
    }
    return(geno_mat[order(geno_mat$Genotype),])
  } else {
    return(data.frame(Genotype = label, Freq = 1))
  }
}

reduce_gcd <- function(num){
  seqs <- purrr::map(as.list(num), ~ seq_len(ceiling(.x/2)) %>%
                c(.x) %>%
                unique)
  divisors <- purrr::map2(seqs, as.list(num), ~ .x[(.y %% .x) == 0])
  gcd <- unlist(divisors) %>%
    table %>%
    as.data.frame() %>%
    dplyr::filter(Freq == length(num)) %>%
    dplyr::select(-Freq) %>%
    unlist() %>%
    as.character() %>%
    as.numeric() %>%
    max()
  num/gcd
}

validate_cross <- function(geno, parent1, parent2, ploidy, reference){
  parents <- furrr::future_map2(parent1, parent2, ~ sort(c(.x, .y)) %>% paste(collapse = "_")) %>%
    unlist()
  ratios <- furrr::future_map(parents, ~ reference[names(reference) == .] %>%
                         magrittr::extract2(1) %>%
                         na.omit)
  furrr::future_map2(geno, ratios, ~ !(.x %in% na.omit(.y)$Genotype)) %>%
    unlist() %>%
    sum()
  # return(ratios)
}

check_crosses <- function(geno, breedr, ploidy=NULL, nc=1){
  future::plan(multisession, workers = nc)
  genos <- colnames(geno)[-1:-3]
  in_breedr <- dplyr::filter(breedr, Ind %in% genos[genos %in% breedr$Ind]) %>%
    na.omit() %>%
    dplyr::filter(Par1 %in% genos) %>%
    dplyr::filter(Par2 %in% genos)
  print(paste0("Complete genotypic crossing data available for: ", stringr::str_flatten(in_breedr$Ind, ", ")))
  print("Running validation of crosses...")
  reference <- geno_table(ploidy)
  anomolies <- c()
  for(i in 1:nrow(in_breedr)){
    cross_geno <- data.frame(Ind = dplyr::select(geno, in_breedr$Ind[i]),
                             Par1 = dplyr::select(geno, in_breedr$Par1[i]),
                             Par2 = dplyr::select(geno, in_breedr$Par2[i])) %>%
      na.omit()
    anomoly <- data.frame(Ind = in_breedr[i,1], Par1 = in_breedr[i,2], Par2 = in_breedr[i,3],
                          Anomolous_Loci = validate_cross(cross_geno[,1], cross_geno[,2], cross_geno[,3], 4, reference))
    anomolies <- rbind(anomolies, anomoly)
    print(paste0(anomolies[i,4],
                 " anomolous loci detected in ", in_breedr$Ind[i]))
  }
  return(anomolies)
}

check_list <- function(test_list, geno, breedr, ploidy, nc=1){
  future::plan(multisession, workers = nc)
  genos <- colnames(geno)[-1:-3]
  print("Running validation of crosses...")
  reference <- geno_table(ploidy)
  anomolies <- c()
  for(i in 1:nrow(test_list)){
    cross_geno <- data.frame(Ind = dplyr::select(geno, test_list$Ind[i]),
                             Par1 = dplyr::select(geno, test_list$Par1[i]),
                             Par2 = dplyr::select(geno, test_list$Par2[i])) %>%
      na.omit()
    anomoly <- data.frame(Ind = test_list[i,1], Par1 = test_list[i,2], Par2 = test_list[i,3],
                          Anomolous_Loci = validate_cross(cross_geno[,1], cross_geno[,2], cross_geno[,3], 4, reference))
    anomolies <- rbind(anomolies, anomoly)
    print(paste0(anomolies[i,4],
                 " anomolous loci detected in ", test_list$Ind[i],
                 " with parents ", test_list$Par1[i], " and ", test_list$Par2[i]))
  }
  return(anomolies)
}

geno_table <- function(ploidy){
  parent_genos <- combn(0:ploidy, ploidy/2) %>%
    cbind(rbind(0:ploidy, 0:ploidy)) %>%
    as.data.frame()
  seg_ratios <- purrr::map(parent_genos, ~ seg_ratio(.[1], .[2], ploidy))
  names(seg_ratios) <- purrr::map(parent_genos, ~ stringr::str_flatten(., collapse = "_")) %>%
    unlist()
  return(seg_ratios)
}

#' Computes heterozygosity
#'
#' This function takes a vector of SNP calls at and returns the heterozygosity of that individual
#'
#' @param geno_vec A vector of either integer or character values indicating the form
#' of the allele
#' @param ploidy A positive integer indicating the ploidy level of the parents and progeny.  Only required if
#' parental genotypes are integers.
#' @return A numeric indicating the proportion of loci that are heterozygous
#' @export

heterozygosity <- function(geno_vec, ploidy){
  geno_vec = na.omit(geno_vec)
  if(is.numeric(geno_vec)){
    H = sum(geno_vec != 0 & geno_vec != ploidy)/length(geno_vec)
  } else if(is.character(geno_vec)){
    H = (purrr::map_int(stringr::str_split(geno_vec, ""), ~ unique(.) %>%
                          length) > 1) %>%
      sum() %>%
      `/` (length(geno_vec))
  }
  return(H)
}

#' Sampling random progeny of parents
#'
#' This function takes a pair of parent SNPs and outputs a randomly sampled offspring
#'
#' @param p1 A numeric or character indicating the first parental genotype
#' @param p2 A numeric or character indicating the second parental genotype
#' @param ploidy A positive integer indicating the ploidy level of the parents and progeny.
#' @return a randomly sampled progeny genotype
#' @export

possible_geno <- function(p1, p2, ploidy){
  geno_prob <- seg_ratio(p1, p2, ploidy) %>%
    mutate(Freq = Freq / sum(Freq, na.rm = T)) %>%
    na.omit()
  if(nrow(geno_prob) > 1){
    o1 <- sample(geno_prob$Genotype, 1, F, prob = geno_prob$Freq)
  } else {
    o1 <- geno_prob$Genotype
  }
  return(o1)
}



###################################################################################################
#  Functions applying segregation ratios to multiple loci
###################################################################################################

multilocus <- function(ploidy, loci){
  alleles1 <- purrr::map(letters[1:loci], ~ str_flatten(rep(., ploidy/2)))
  alleles2 <- purrr::map(LETTERS[1:loci], ~ str_flatten(rep(., ploidy/2)))
  parents <- paste0(alleles1, alleles2)
  seg <- purrr::map(parents, ~ seg_ratio(.)$Genotype) %>%
    purrr::cross() %>%
    purrr::map(~ stringr::str_flatten(.)) %>%
    unlist()
  return(seg)
}

#Precalculating an index for genotype input options
# multilocus_table <- list()
# for (i in 1:4) {
#   multilocus_table[[i]] <- map(c(2, 4, 6, 8, 10), ~ multilocus(., i))
#   names(multilocus_table[[i]]) <- c(2, 4, 6, 8, 10)
# }
# use_data(multilocus_table, internal = T, overwrite = T)

#' Expected segregation ratios considering multiple loci.
#'
#' This currently a low-level function called by xplotter that is not intended for direct use
#'
#' @export

multi_seg <- function(parent1, parent2, ploidy, loci){
  parent1_sep <- substring(parent1, seq(1, ((ploidy*loci)-(ploidy-1)), ploidy),
                       seq(ploidy, ploidy*loci, ploidy))
  parent2_sep <- substring(parent2, seq(1, ((ploidy*loci)-(ploidy-1)), ploidy),
                       seq(ploidy, ploidy*loci, ploidy))
  ratios <- purrr::map2(parent1_sep, parent2_sep, ~ seg_ratio(.x, .y) %>%
                   na.omit %>%
                   dplyr::mutate(Freq = (Freq/sum(Freq))))
  classes <- ratios[[1]]$Genotype
  freqs <- ratios[[1]]$Freq
  if(loci > 1){
    for(i in 2:(length(ratios))){
      classes <- purrr::cross2(classes, ratios[[i]]$Genotype) %>%
        purrr::map(~ stringr::str_flatten(.)) %>%
        unlist()
      freqs <- purrr::cross2(freqs, ratios[[i]]$Freq) %>%
        purrr::map(~ prod(unlist(.))) %>%
        unlist()
    }
  }
  return(data.frame(Genotype = classes, Freq = freqs))
}
