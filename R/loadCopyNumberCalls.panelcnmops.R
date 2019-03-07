#' loadCopyNumberCalls.panelcnmops
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
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls.panelcnmops(panelcn.mops.res, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = NA, segment.value.col = NULL, genome = NULL, verbose = TRUE)
#' 
#' @param panelcn.mops.res The name of the file with the data or the name of the variable with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to NULL)
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with CN information. If NULL, it is automatically identified. (default to NA)
#' @param segment.value.col The name or number of the column with segment value. If NULL, it is automatically identified. (default to NULL)
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#' ## loadCopyNumberCalls.panelcnmops from panelcn.mops resulttable:
#' data(panelcn.mops, package = "panelcn.mops")
#' XandCB <- test
#' sampleNames <- colnames(elementMetadata(XandCB))
#' elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(control))
#' 
#' resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, countWindows = countWindows, sampleNames = sampleNames)
#' panelcn.mops.res <- resulttable[[1]]
#' loadCopyNumberCalls.panelcnmops(panelcn.mops.res = panelcn.mops.res)
#' 
#' @export loadCopyNumberCalls.panelcnmops
#'

loadCopyNumberCalls.panelcnmops <- function(panelcn.mops.res, 
                                            chr.col = NULL, 
                                            start.col = NULL,
                                            end.col = NULL, 
                                            cn.col = NA,
                                            segment.value.col = NULL, 
                                            genome = NULL, 
                                            verbose=TRUE){

 
  segs <- loadCopyNumberCalls(cnv.data = panelcn.mops.res, 
                              chr.col= chr.col,  
                              start.col = start.col,  
                              end.col = end.col,  
                              cn.col = cn.col,  
                              segment.value.col = segment.value.col, 
                              genome = genome, 
                              verbose = verbose) 
  
  GenomicRanges::mcols(segs)$cn <- as.integer(gsub(x = GenomicRanges::mcols(segs)$CN, pattern = "CN", replacement = ""))
  
  GenomeInfoDb::seqlevelsStyle(segs) <- "UCSC"
  
  #if there is more than one sample in cn.mops.res
  if(length(table(segs$Sample)) > 1){
    segs <- split(segs, GenomicRanges::mcols(segs)$Sample)
  }
  
  
  return(segs)
                                            
}