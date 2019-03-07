#' loadCopyNumberCalls.CNVkit
#'
#' @description
#' Loads copy number calls from CNVkit.cns file format
#'
#' @details
#' This function will load segments data from CNVkit.cns file format. 
#' Internally it uses the \code{toGRanges} function from regioneR
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
#' @usage loadCopyNumberCalls.CNVkit(cnvkit.file, chr.col = "chromosome", start.col = "start", end.col = "end", cn.col = NA, segment.value.col = "log2", genome = NULL, verbose=TRUE)
#'
#' @param cnvkit.file The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to "chromosome")
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to "start")
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to "end")
#' @param cn.col The name or number of the column with CN information. If NULL, it is automatically identified. (default to NA)
#' @param segment.value.col The name or number of the column with segment value. If NULL, it is automatically identified. (default to "log2")
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#' ## loadCopyNumberCalls.CNVkit from .cns file format:
#' ## An example of .cns file format is found in https://github.com/etal/cnvkit/blob/master/test/formats/cl_seq.cns.
#' cnvkit.file <- "testFiles/cl_seq.cns"
#' loadCopyNumberCalls.CNVkit(cnvkit.file)
#'
#' @export loadCopyNumberCalls.CNVkit
#'

loadCopyNumberCalls.CNVkit <- function(cnvkit.file, 
                                       chr.col = "chromosome", 
                                       start.col = "start",
                                       end.col = "end", 
                                       segment.value.col = "log2", 
                                       cn.col = NA,
                                       genome = NULL,
                                       verbose = TRUE){
  #Load file.seg. First, probe if file exist
  if(!file.exists(cnvkit.file)){
    stop(paste0(cnvkit.file, " does not exist or you are not in the correct directory."))
  }
  
  segs <- loadCopyNumberCalls(cnv.data = cnvkit.file, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NA, 
                              genome = genome,
                              verbose = verbose)
  

  return(segs)
}


