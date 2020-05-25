#' loadCopyNumberCallsCNVkit
#'
#' @description
#' Loads copy number calls from CNVkit.cns file format
#'
#' @details
#' This function will load segments data from CNVkit.cns file format and from CNVkit.cnr. 
#' Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCallsCNVkit(cnvkit.file, chr.col = "chromosome", start.col = "start", end.col = "end", segment.value.col = "log2", file.type = "auto", cn.col = NULL, zero.based = TRUE, genome = NULL, verbose = TRUE)
#' 
#' @param cnvkit.file The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to "chromosome")
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to "start")
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to "end")
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to "log2")
#' @param file.type (character) wheter to load ".cns", ".cnr" or "auto" if the file type is automatically recognised. (defaults to "auto)   
#' @param zero.based (logical) Whether the data is zero-based and half open (i.e. ranges are defined by (start:end] so chr1:10-20 represents nine bases long features spanning from base 11 to 20). (defaults to FALSE)
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
#' # Loading a CNVkit.cns file 
#' cnvkit.file <- system.file("extdata", "CNVkit_output.cns", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCallsCNVkit(cnvkit.file)
#'
#' # Loading a CNVkit.cnr file
#' cnvkit.file <- system.file("extdata", "CNVkit_output.cnr", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCallsCNVkit(cnvkit.file)
#' 
#' @export loadCopyNumberCallsCNVkit
#'

loadCopyNumberCallsCNVkit <- function(cnvkit.file, 
                                      chr.col = "chromosome", 
                                      start.col = "start",
                                      end.col = "end", 
                                      segment.value.col = "log2", 
                                      file.type = "auto",
                                      cn.col = NULL,
                                      zero.based = TRUE,
                                      genome = NULL,
                                      verbose = TRUE){
 
   #Load .cns or .cnr  file. First, probe if file exist
  if(is.character(cnvkit.file)){
    if(!file.exists(cnvkit.file)){
    stop(paste0(cnvkit.file, " does not exist or you are not in the correct directory."))
    }
  }
  
  # test if file.type parameter is a character
  if(is.character(file.type)){
    if(file.type != "auto" && file.type != ".cnr"&& file.type != ".cns"){
    stop("file.type paremeter must be auto, .cnr or .cns")
    }
    
  }else{
    stop("file.type paremeter must be a character. Either auto, .cnr or .cns")
  }

  
  
  # we detect the file type to load
  if(file.type == "auto" ||file.type == ".cnr"||file.type == ".cns"){
    
    if(grepl(pattern = ".cnr", x = cnvkit.file)){
      # segment.value is lrr in .cnr files
      chng.segment.col.name <- "lrr"
    }
    
    else if (grepl(pattern = ".cns", x = cnvkit.file)){
      # segment.value is segement.value in .cns files
      chng.segment.col.name <- "segement.value"
    
    }else{
      stop(paste0(cnvkit.file, " is not a CNVkit.cnr or CNVkit.cns file."))
      
    }
  }
  
  
  segs <- loadCopyNumberCalls(cnv.data = cnvkit.file, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NULL, 
                              genome = genome,
                              zero.based = zero.based,
                              verbose = verbose)
  
  names(mcols(segs))[names(mcols(segs)) == "segment.value"]<- chng.segment.col.name
  
  
  return(segs)
}