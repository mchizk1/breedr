#' Assembling Pedigrees.
#'
#' This function uses UA crossing records to quickly assemble a pedigree matrix.
#'
#' @param breedr A data.frame containing three columns: individual, female parent,
#' and male parent
#' @param genotype A character string indicating the name of the cultivar or breeding
#' selection that you wish to assemble a pedigree for.  This function will attempt
#' to recognize all variations of the name by calling nameR().
#' @param n (default = 4) An positive integer indicating the number of generations
#' you wish to see represented in the pedigree.
#' @param orientation (default = TB) A character string.  Options: "TB" for Top to Bottom,
#' "LR" for Left to Right, "RL" for Right to Left
#' @param method (default = FULL) A string indicating the method for plotting a
#' pedigree.  Options: 'FULL' for a full pedigree with parental sex color formatting.
#' 'CA' for grouped Common Ancestors.
#' @param str_ops A character vector containing optional methods for string manipulation.
#' Possible options include: \cr
#' 'strip_ws' (default) - Strip leading and trailing white space. \cr
#' 'upper' (default) - Make all uppercase to eliminate case differences.\cr
#' 'lower' - Make all lowercase to eliminate case differences.  Cannot be combined with 'upper'. \cr
#' 'shrink_ws' (default) - remove duplicated white space from center of strings. \cr
#' NULL is acceptable to skip string manipulation
#' @param na_val An optional string indicating a value to be treated as missing.
#' @param color1 The primary color. In the "FULL" method, this is the male parent
#' and in the "CA" method this is the common ancestors.
#' @param color2 The secondary color. In the "FULL" method, this is the female parent
#' and in the "CA" method, this is the non-common ancestors.
#' @param color3 The tertiary color. In both methods, this is the color of the individual
#' indicated in the \code{genotype} argument
#' @return A plot of assembled pedigree going back \code{n} generations.
#' @export
#' @examples
#' plotigree(habsburg, "Ferdinand III HRE")
#' plotigree(habsburg, "Charles the Bewitched", method = "CA")

plotigree <- function(breedr, genotype, n = 4, orientation = "TB", method = "FULL",
                      str_ops = c("strip_ws", "upper", "shrink_ws"), na_val = NULL,
                      color1 = NULL, color2 = NULL, color3 = NULL){
  breedr <- new_breedr(breedr, str_ops, na_val)
  assertthat::assert_that(assertthat::is.string(genotype),
                          msg = "genotype must be a character string")
  assertthat::assert_that(orientation %in% c("TB", "LR", "RL"),
                          msg = "orientation must be one of the following: TB, LR, or RL")
  genotype <- nameR(genotype, str_ops, na_val)
  assertthat::assert_that(genotype %in% breedr$Ind,
                          msg = paste0(genotype, " was not found in the provided breedr object"))
  assertthat::assert_that(assertthat::is.count(n),
                          msg = "n must be a positive integer")
  assertthat::assert_that(method %in% c("FULL", "CA"),
                          msg = "method must be either 'FULL' or 'CA'")
  if (method == "FULL"){
    dot_script <- FULL_method(breedr, genotype, n, orientation, color1 = color1,
                              color2 = color2, color3 = color3)
  } else if (method == "CA"){
    dot_script <- CA_method(breedr, genotype, n, orientation, color1 = color1,
                            color2 = color2, color3 = color3)
  }
  DiagrammeR::grViz(dot_script)
  #return(dot_script)
}

#===============================================================================
# This function takes information about a node on a dot flowchart and writes
# three lines of dot code defining the node, node label, and edge connections.
#===============================================================================

writedot <- function(nodeid, label=NULL, attachid, color){
  if(!is.numeric(nodeid)){
    nodeid <- paste0("'",nodeid,"'")
  }
  if(!is.numeric(attachid)){
    attachid <- paste0("'",attachid,"'")
  }
  if(!is.null(label)){
    # Assembling DOT script node IDs
    nodedef <- paste0(
      (nodeid)," [label = '@@",(nodeid),"', style = filled, color = ",color,"] \n"
    )
    # Assembling DOT script node labels
    nodelab <- paste0(
      "[",nodeid,"]: '",label,"' \n"
    )
    # Assembling DOT script edge definitions
    edgedef <- paste0(
      (nodeid)," -> ",attachid, "\n")
    dotvec <- c(nodedef, nodelab, edgedef)
  } else {
    # Assembling DOT script node IDs
    nodedef <- paste0(
      (nodeid)," [style = filled, color = ",color,"] \n"
    )
    # Assembling DOT script edge definitions
    edgedef <- paste0(
      (nodeid)," -> ",attachid, "\n")
    dotvec <- c(nodedef, "", edgedef)
  }
  return(dotvec)
}

#===============================================================================
# This function writes a dot header and assembles all of the components of the
# final dot script for grViz.
#===============================================================================

assembledot <- function(nodedef, nodelab=NULL, edgedef, orientation){
  if(!is.null(nodelab)){
    dot_start <- paste0("\n digraph pedigree{ \n
    graph [rankdir = ",orientation,"] \n
    node [shape = box, fontname = Helvetica] \n")
    fulldot <- paste0(dot_start, nodedef,
                      edgedef, "} \n",
                      nodelab)
  } else {
    dot_start <- paste0("\n digraph pedigree{ \n
    graph [rankdir = ",orientation,"] \n
    node [shape = box, fontname = Helvetica] \n")
    fulldot <- paste0(dot_start, nodedef,
                      edgedef, "}")
  }
  return(fulldot)
}

#===============================================================================
# This function supports the CA method of pedigree plotting
#===============================================================================

CA_method <- function(breedr, genotype, n, orientation = "TB", color1 = color1,
                      color2 = color2, color3 = color3){
  if(is.null(color1)){
    color1 = "Khaki"
  } else {
    assertthat::assert_that(assertthat::is.string(color1), msg = "color1 must be a string")
  }
  if(is.null(color2)){
    color2 = "Gray"
  } else {
    assertthat::assert_that(assertthat::is.string(color2), msg = "color2 must be a string")
  }
  if(is.null(color3)){
    color3 = "Khaki"
  } else {
    assertthat::assert_that(assertthat::is.string(color2), msg = "color3 must be a string")
  }
  # Assembling data needed for CA method of plotting
  parents <- dplyr::filter(breedr, Ind == genotype)
  P1 <- ped_idx(parents$Par1, breedr)
  P2 <- ped_idx(parents$Par2, breedr)
  P1_anc <- unique(P1$pedvec)
  P2_anc <- unique(P2$pedvec)
  ca <- ca_detect(P1_anc, P2_anc)
  Pboth <- ped_idx(genotype, breedr)
  if(!is.null(ca)){
    # Shrink vectors indices based on shared ancestry
    for(i in ca){
      Pboth <- distinct_ca(Pboth, i, breedr)
      P1 <- distinct_ca(P1, i)
      P2 <- distinct_ca(P2, i)
    }
  }
  ca <- unique(ca_detect(P1$pedvec, P2$pedvec))
  # Now it's time to use the assembled data to write a DOT script
  nodedef <- paste0("'",genotype,"' [style=filled, color=",color3,"] \n")
  edgedef <- c("\n")
  defined <- c()
  # What happens if the individual resulted from a self?
  if (parents[1,2] == parents[1,3]){
    dotvec <- writedot(nodeid = parents[1,2], attachid = genotype, color = color1)
    nodedef <- paste0(nodedef, dotvec[1])
    edgedef <- paste0(edgedef, dotvec[3])
  } else {
    # Otherwise, for each generation in the parent indices...
    #for (i in 1:max(Pboth$gidx)){
    for (i in 1:n){
      active_gen <- dplyr::filter(Pboth, gidx == i)
      # And for each occupied position in that generation...
      for (j in active_gen$pidx){
        active_pos <- dplyr::filter(active_gen, pidx == j)
        # different edge definitions based on odd or even position index
        if((j %% 2) == 0){
          attach_pos <- j/2
        } else {
          attach_pos <- (j+1)/2
        }
        attachid <- dplyr::filter(Pboth, gidx == i-1 & pidx == attach_pos)
        if(active_pos$pedvec %in% ca){
          color <- color1
        } else {
          color <- color2
        }
        if(!(active_pos$pedvec %in% defined)){
          dotvec <- writedot(nodeid = active_pos$pedvec, attachid = attachid, color = color)
          nodedef <- paste0(nodedef, dotvec[1])
          if (!grepl(dotvec[3], edgedef)){
            edgedef <- paste0(edgedef, dotvec[3])
          }
          defined <- paste0(defined, active_pos$pedvec)
        } else {
          dotvec <- writedot(nodeid = active_pos$pedvec, attachid = attachid, color = color)
          if (!grepl(dotvec[3], edgedef)){
            edgedef <- paste0(edgedef, dotvec[3])
          }
        }
      }
    }
  }
  dotscript <- assembledot(nodedef = nodedef, edgedef = edgedef,
                           orientation = orientation)
  return(dotscript)
}


#===============================================================================
# This function supports the FULL method of pedigree plotting
#===============================================================================

FULL_method <- function(breedr, genotype, n, orientation = "TB", color1 = color1,
                        color2 = color2, color3 = color3){
  if(is.null(color1)){
    color1 = "PaleTurquoise"
  } else {
    assertthat::assert_that(assertthat::is.string(color1), msg = "color1 must be a string")
  }
  if(is.null(color2)){
    color2 = "Pink"
  } else {
    assertthat::assert_that(assertthat::is.string(color2), msg = "color2 must be a string")
  }
  if(is.null(color3)){
    color3 = "Thistle"
  } else {
    assertthat::assert_that(assertthat::is.string(color2), msg = "color3 must be a string")
  }
  nodecount <- 2
  nodebase <- 1
  nodeloop <- 1
  attachid <- 1
  nodeloc <- c()
  pvec_old <- c(genotype)
  pvec_new <- c()
  label_idx <- data.frame(label=c(), node=c())
  gensz <-1
  nodedef<-paste0((1)," [label = '@@1', style=filled, color=",color3," ] \n")
  nodelab<-paste0("[",1,"]: '",genotype,"' \n")
  edgedef<-c()
  # This loop was created to trace back pedigrees using the requested number of
  # generations
  # for each generation requested
  for(i in 1:n){
    k <- 1
    # for each individual that had parental data
    for(j in nodeloop){
      parents <- dplyr::filter(breedr, Ind == pvec_old[j])
      colnames(parents) <- NULL
      pvec_new <- c(pvec_new, unlist(parents[1,2:3]))
      if(!is.na(parents[1,2]) & !is.na(parents[1,3])){
        dotvec <- writedot(nodecount, parents[1,2], attachid[j], color = color2)
        nodedef <- paste0(nodedef, dotvec[1])
        nodelab <- paste0(nodelab, dotvec[2])
        edgedef <- paste0(edgedef, dotvec[3])
        nodecount=nodecount+1
        dotvec <- writedot(nodecount, parents[1,3], attachid[j], color = color1)
        nodedef <- paste0(nodedef, dotvec[1])
        nodelab <- paste0(nodelab, dotvec[2])
        edgedef <- paste0(edgedef, dotvec[3])
        nodecount=nodecount+1
        nodeloc <- c(nodeloc, k, k+1)
        pvec <- c(parents[1,2:3])
        k=k+2
      } else if (!is.na(parents[1,2]) & is.na(parents[1,3])){
        dotvec <- writedot(nodecount, parents[1,2], attachid[j], color = color2)
        nodedef <- paste0(nodedef, dotvec[1])
        nodelab <- paste0(nodelab, dotvec[2])
        edgedef <- paste0(edgedef, dotvec[3])
        nodecount=nodecount+1
        nodeloc <- c(nodeloc, k)
        k=k+1
      } else if (is.na(parents[1,2]) & !is.na(parents[1,3])){
        dotvec <- writedot(nodecount, parents[1,3], attachid[j], color = color1)
        nodedef <- paste0(nodedef, dotvec[1])
        nodelab <- paste0(nodelab, dotvec[2])
        edgedef <- paste0(edgedef, dotvec[3])
        nodecount=nodecount+1
        nodeloc <- c(nodeloc, k)
        k=k+1
      } else {
      }
    }
    nodeloop <- nodeloc
    attachid <- nodeloop+nodebase
    nodebase <- nodecount-1
    gensz=gensz*2
    nodeloc <- c()
    pvec_old <- pvec_new[!is.na(pvec_new)]
    pvec_new <- c()
  }
  dot_script <- assembledot(nodedef, nodelab, edgedef, orientation)
  return(dot_script)
}
