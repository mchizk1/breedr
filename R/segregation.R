
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
    sort()
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
            sort %>%
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

heterozygosity <- function(geno_vec){
  no_na <- na.omit(geno_vec)
  (!(no_na %in% c(0,4))) %>%
    sum() %>%
    `/`(length(no_na))
}

