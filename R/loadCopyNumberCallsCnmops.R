#' loadCopyNumberCallsCnmops
#'
#' @description
#' Loads copy number calls from either cn.mops result class or in a tabular format
#'
#' @details
#' This function will load segments data from either cn.mops result class or
#' in a tabular format from a file. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' Column names of the format loaded are specified as default but,
#' if no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCallsCnmops(cn.mops.res, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = "CN", segment.value.col = "median", genome = NULL, verbose = TRUE)
#'
#' @param cn.mops.res The name of the file with the data or the name of the variable with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to "CN")
#' @param segment.value.col (number or character) The name or number of the column with segment value. It can be mean or median. If NULL, it is automatically identified. (defaults to "median")
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment or a list of GRanges with a GRanges per sample.
#'
#' @examples
#' # loadCopyNumberCallsCnmops from cn.mops result class:
#' #NOT RUN - Example loading data from a results object created by cn.mops
#' # require(cn.mops)
#' # data(cn.mops, package = "cn.mops")
#' 
#' # cn.mops.res <- cn.mops(XRanges)
#' # cn.mops.res <- calcIntegerCopyNumbers(cn.mops.res)
#' # cnv.call <- loadCopyNumberCallsCnmops(cn.mops.res = cn.mops.res)  
#' 
#' # loadCopyNumberCallsCnmops from a file where cn.mops result was saved:
#' cn.mops.res <- system.file("extdata", "cn.mops.segmentation.csv", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCallsCnmops(cn.mops.res = cn.mops.res) 
#'
#' @export loadCopyNumberCallsCnmops
#'
#' @importFrom cn.mops segmentation
#' 

loadCopyNumberCallsCnmops<- function(cn.mops.res, 
                                      chr.col = NULL, 
                                      start.col = NULL,
                                      end.col = NULL, 
                                      cn.col = "CN",
                                      segment.value.col = "median", 
                                      genome = NULL, 
                                      verbose=TRUE){
  
  #if cn.mops.res is "CNVDetectionResult" we need to extract the results using segmentation from cn.mops package
  if(methods::is(cn.mops.res, "CNVDetectionResult")){
    cn.mops.res <- cn.mops::segmentation(cn.mops.res)
  }  
  
  segs <- loadCopyNumberCalls(cnv.data = cn.mops.res, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col, 
                              loh.col = NULL,
                              genome = genome,
                              verbose = verbose)
  
  GenomicRanges::mcols(segs)$cn <- as.integer(gsub(segs$cn, pattern = "CN", replacement = ""))
  
  if(any(c("width", "strand") %in% names(GenomicRanges::mcols(segs)))){
    GenomicRanges::mcols(segs) <-  GenomicRanges::mcols(segs)[,!(names(GenomicRanges::mcols(segs)) %in% c("width", "strand"))]
    
  }
  
  
  #if there is more than one sample in cn.mops.res
  if(length(unique(segs$sampleName)) > 1){
    segs <- split(segs, mcols(segs)$sampleName)
  }
  
  return(segs)
}  
