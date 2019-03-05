#' loadCopyNumberCalls.DECoN
#'
#' @description
#' Loads copy number calls from a tabular format
#'
#' @details
#' This function will load segments data from any "bed-like" data structure
#' in R or file. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files.
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls(cnv.data, chr.col=NULL, start.col=NULL, end.col=NULL, cn.col=NULL, loh.col=NULL, segment.value.col=NULL, genome=NULL, verbose=TRUE)
#'
#' @param cnv.data The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to NULL)
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with BAF information. If NULL, it is automatically identified. (default to NULL)
#' @param loh.col The name or number of the column with LRR information. If NULL, it is automatically identified. (default to NULL)
#' @param segment.value.col The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (default to NULL)
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#'
#'
#'
#' @export loadCopyNumberCalls.DECoN
#'
#' @importFrom GenomicRanges mcols
#' @importFrom GenomeInfoDb seqlevelsStyle


cnv.data <- decon.file
segs <- 0
#loadCopyNumberCalls.DECoN

loadCopyNumberCalls.DECoN <- function(decon.file, 
                                      chr.col = NULL, 
                                      start.col = "Start",
                                      end.col = "End",
                                      cn.col = NA,
                                      segment.value.col= NULL,
                                      genome = NULL, 
                                      verbose = TRUE){
  #load data
  #Load decon.file First, probe if file exist
  if(!file.exists(decon.file)){
    stop(paste0(decon.file, " does not exist or you are not in the correct directory."))
  }
  
  segs <- loadCopyNumberCalls(cnv.data = decon.file,
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NA, 
                              genome = genome,
                              verbose = verbose)
  # We order the segments
  
  
  # Adding a new column (cn) where CNV.type information is translated in an integer
  # depending on whether there is a duplication (3) or a deletion (1)
  GenomicRanges::mcols(segs)$cn[GenomicRanges::mcols(segs)$CNV.type == "duplication"]<- 3
  GenomicRanges::mcols(segs)$cn[GenomicRanges::mcols(segs)$CNV.type == "deletion"]<- 1
  
  #if there is more than one sample in cn.mops.res
  if(length(table(segs$Sample)) > 1){
    segs <- split(segs, mcols(segs)$Sample)
    segs <- lapply(segs, function(segments.to.sort) sort(segments.to.sort))
    
  }else{
    segs <- sort(segs)
  }

  return(segs)
}

#Example:

decon.file <- "testFiles/calls_all.txt"
cc<-loadCopyNumberCalls.DECoN(decon.file = decon.file, genome = "hg19")
  