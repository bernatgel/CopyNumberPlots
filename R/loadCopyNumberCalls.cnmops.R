#' loadCopyNumberCalls.cnmops
#'
#' @description
#' Loads copy number calls from either cn.mops result class or in a tabular format
#'
#' @details
#' This function will load segments data from either cn.mops result class or
#' in a tabular format or file. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' Column names of the format loaded are specified as default but,
#' if no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls.cnmops(cn.mops.res, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = "CN", segment.value.col = "median", genome = NULL, verbose = TRUE)
#'
#' @param cn.mops.res The name of the file with the data or the name of the variable with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to NULL)
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with CN information. If NULL, it is automatically identified. (default to "CN")
#' @param segment.value.col The name or number of the column with segment value. It can be mean or median. If NULL, it is automatically identified. (default to "median")
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#' ## loadCopyNumberCalls.cnmops from cn.mops result class:
#' data(cn.mops, package = "cn.mops")
#' 
#' cn.mops.res <- cn.mops(XRanges)
#' cn.mops.res <- calcIntegerCopyNumbers(cn.mops.res)
#' loadCopyNumberCalls.cnmops(cn.mops.res = cn.mops.res)  
#' 
#' ## loadCopyNumberCalls.cnmops from a file where cn.mops result was saved:
#' cn.mops.res <- "cn.mops.segmentation.csv"
#' loadCopyNumberCalls.cnmops(cn.mops.res = cn.mops.res) 
#'
#' @export loadCopyNumberCalls.cnmops
#'
#' @importFrom cn.mops segmentation

loadCopyNumberCalls.cnmops<- function(cn.mops.res, 
                                      chr.col = NULL, 
                                      start.col = NULL,
                                      end.col = NULL, 
                                      cn.col = "CN",
                                      segment.value.col = "median", 
                                      genome = NULL, 
                                      verbose=TRUE){
  
  #if cn.mops.res is "CNVDetectionResult"
  if(methods::is(cn.mops.res, "CNVDetectionResult")){
    segs <- cn.mops::segmentation(cn.mops.res)
    segs$cn <- as.integer(gsub(segs$CN,pattern = "CN",replacement = ""))
    segs <- loadCopyNumberCalls(cnv.data = segs,
                                chr.col = chr.col, 
                                start.col = start.col, 
                                end.col = end.col, 
                                cn.col = cn.col, 
                                segment.value.col = segment.value.col,
                                loh.col = NA,
                                genome = genome,
                                verbose = verbose)

  }else{
    segs <- loadCopyNumberCalls(cnv.data = cn.mops.res, 
                                chr.col = chr.col, 
                                start.col = start.col,
                                end.col = end.col, 
                                cn.col = cn.col, 
                                segment.value.col = segment.value.col, 
                                loh.col = NA,
                                genome = genome,
                                verbose = verbose)
    
    GenomicRanges::mcols(segs)$cn <- as.integer(gsub(segs$cn,pattern = "CN",replacement = ""))
    GenomicRanges::mcols(segs) <-  GenomicRanges::mcols(segs)[,-which(names(GenomicRanges::mcols(segs)) %in% c("width","strand"))]
  }
  
  #if there is more than one sample in cn.mops.res
  if(length(table(segs$sampleName)) > 1){
     segs <- split(segs, mcols(segs)$sampleName)
  }
  
  return(segs)
}  
