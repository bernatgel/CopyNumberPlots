#' loadCopyNumberCallsCNVkit
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
#' @usage loadCopyNumberCallsCNVkit(cnvkit.file, chr.col = "chromosome", start.col = "start", end.col = "end", segment.value.col = "log2", cn.col = NULL, genome = NULL, verbose = TRUE)
#' 
#' @param cnvkit.file The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to "chromosome")
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to "start")
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to "end")
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to "log2")
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment. 
#'
#' @examples
#' ## loadCopyNumberCallsCNVkit from .cns file format:
#' ## An example of .cns file format is found at https://github.com/etal/cnvkit/blob/master/test/formats/cl_seq.cns.
#' 
#' cnvkit.file <- system.file("extdata", "CNVkit_output.cns", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCallsCNVkit(cnvkit.file)
#'
#' @export loadCopyNumberCallsCNVkit
#'

loadCopyNumberCallsCNVkit <- function(cnvkit.file, 
                                       chr.col = "chromosome", 
                                       start.col = "start",
                                       end.col = "end", 
                                       segment.value.col = "log2", 
                                       cn.col = NULL,
                                       genome = NULL,
                                       verbose = TRUE){
 
   #Load .seg file. First, probe if file exist
  if(is.character(cnvkit.file)){
    if(!file.exists(cnvkit.file)){
    stop(paste0(cnvkit.file, " does not exist or you are not in the correct directory."))
    }
    
    #If is a .cns file this means we will have segments and/or cn data zero based.
    if(grepl(pattern = ".cnr", x = cnvkit.file)) zero.based  <- FALSE
    if(grepl(pattern = ".cns", x = cnvkit.file)) zero.based <- TRUE
  }
  
  
  segs <- loadCopyNumberCalls(cnv.data = cnvkit.file, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NULL, 
                              genome = genome,
                              zero.based = zero.based
                              verbose = verbose)
  
  
  return(segs)
}


