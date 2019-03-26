#' loadCopyNumberCalls.seg
#'
#' @description
#' Loads copy number calls from .seg file format
#'
#' @details
#' This function will load segments data from .seg file format. 
#' Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCalls.seg (seg.file, chr.col = "chrom", start.col = "loc.start", end.col = "loc.end", segment.value.col = "seg.mean", cn.col = NULL, genome = NULL, chr.transformation = "23:X,24:Y,25:MT", verbose = TRUE)
#'
#' @param seg.file The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to "chrom")
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to "loc.start")
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to "loc.end")
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to "seg.mean")
#' @param cn.col (number or character)The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param chr.transformation (character)(character) The transformation of the chromosome names in a comma separated "key:value" format as detailed at https://cnvkit.readthedocs.io/en/stable/importexport.html#import-seg.(defaults to "23:X,24:Y,25:MT")
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment
#'
#' @examples
#' ## loadCopyNumberCalls.seg from .seg file format: 
#' ## the file to run in the example can be found in: https://software.broadinstitute.org/software/igv/SEG
#' ## under example.seg file name.
#' 
#' seg.file <- system.file("extdata", "DNACopy_output.seg", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCalls.seg(seg.file = seg.file)
#'
#' @export loadCopyNumberCalls.seg
#'

loadCopyNumberCalls.seg <- function(seg.file, 
                                    chr.col = "chrom", 
                                    start.col = "loc.start",
                                    end.col = "loc.end", 
                                    segment.value.col = "seg.mean", 
                                    cn.col = NULL,
                                    genome = NULL,
                                    chr.transformation = "23:X,24:Y,25:MT",
                                    verbose = TRUE){
  
  #Load file.seg. First, probe if file exist
  if(!file.exists(seg.file)) stop(paste0(seg.file, " does not exist or you are not in the correct directory."))
  
  header <- strsplit(x = gsub(readLines(seg.file, n = 1), pattern = "'", replacement = ""), split = "\t")[[1]]
  
  seg.data <- read.table(file = seg.file, sep = "\t", skip = 1, stringsAsFactors = FALSE)
  colnames(seg.data) <-  header
  
  if(!is.null(chr.transformation)){
    chr.col <- getChrColumn(df = seg.data, col = chr.col, verbose = FALSE)
    seg.data[,chr.col] <- transformChr(chr = seg.data[,chr.col], chr.transformation = chr.transformation)
  }
  
  
  segs <- loadCopyNumberCalls(cnv.data = seg.data, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NULL, 
                              genome = genome,
                              verbose = verbose)
  

  #if there is more than one sample in cn.mops.res
  if(length(table(segs$ID)) > 1){
    segs <- split(segs, mcols(segs)$ID)
  }
  
  return(segs)
}
