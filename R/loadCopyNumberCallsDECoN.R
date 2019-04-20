#' loadCopyNumberCallsDECoN
#'
#' @description
#' Loads copy number calls from DECoN  output file
#'
#' @details
#' This function will load segments data from DECoN output file. 
#' Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCallsDECoN(decon.file, chr.col = NULL, start.col = "Start", end.col = "End", cn.col = NULL, segment.value.col = NULL, genome = NULL, verbose = TRUE)
#'
#' @param decon.file The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to "Start")
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to "End")
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NA)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to NULL)
#' @param genome (character) The name of the genome. (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment or a list of GRanges with a GRanges per sample.
#'
#' @examples
#' decon.file <- system.file("extdata", "DECoN_output.txt", package = "CopyNumberPlots", mustWork = TRUE)
#' cn.calls <- loadCopyNumberCallsDECoN(decon.file = decon.file)
#'
#' @export loadCopyNumberCallsDECoN
#'
#' @importFrom GenomicRanges mcols
#' @importFrom GenomeInfoDb seqlevelsStyle

loadCopyNumberCallsDECoN <- function(decon.file, 
                                      chr.col = NULL, 
                                      start.col = "Start",
                                      end.col = "End",
                                      cn.col = NULL,
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
                              loh.col = NULL, 
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
