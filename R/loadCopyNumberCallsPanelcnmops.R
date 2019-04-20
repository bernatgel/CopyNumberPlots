#' loadCopyNumberCallsPanelcnmops
#' 
#' @description
#' Loads copy number calls from either cn.mops result class or in a tabular format
#'
#' @details
#' This function will load segments data from panelcn.mops resulttable or
#' in a tabular format or file. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCallsPanelcnmops(panelcn.mops.res, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = NULL, segment.value.col = NULL, genome = NULL, verbose = TRUE)
#' 
#' @param panelcn.mops.res The name of the file with the data or the name of the variable with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to NULL)
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (character) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment or a list of GRanges with a GRanges per sample.
#'
#' @examples
#' ## loadCopyNumberCallsPanelcnmops from panelcn.mops resulttable:
#' library(panelcn.mops)
#' data(panelcn.mops, package = "panelcn.mops")
#' XandCB <- test
#' sampleNames <- colnames(elementMetadata(XandCB))
#' elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(control))
#' 
#' resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, countWindows = countWindows, sampleNames = sampleNames)
#' panelcn.mops.res <- resulttable[[1]]
#' cnv.call <- loadCopyNumberCallsPanelcnmops(panelcn.mops.res = panelcn.mops.res)
#' 
#' @export loadCopyNumberCallsPanelcnmops
#'

loadCopyNumberCallsPanelcnmops <- function(panelcn.mops.res, 
                                            chr.col = NULL, 
                                            start.col = NULL,
                                            end.col = NULL, 
                                            cn.col = NULL,
                                            segment.value.col = NULL, 
                                            genome = NULL, 
                                            verbose = TRUE){

 
  segs <- loadCopyNumberCalls(cnv.data = panelcn.mops.res, 
                              chr.col= chr.col,  
                              start.col = start.col,  
                              end.col = end.col,  
                              cn.col = cn.col,  
                              segment.value.col = segment.value.col, 
                              genome = genome, 
                              verbose = verbose) 

  GenomicRanges::mcols(segs)$cn <- as.integer(gsub(x = GenomicRanges::mcols(segs)$cn, pattern = "CN", replacement = ""))
  
  
  #if there is more than one sample in cn.mops.res
  if(length(unique(segs$Sample)) > 1){
    segs <- split(segs, GenomicRanges::mcols(segs)$Sample)
  }
  
  
  return(segs)
  
}