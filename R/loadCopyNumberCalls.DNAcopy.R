#' loadCopyNumberCalls.DNAcopy
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
#' @usage loadCopyNumberCalls.DNAcopy(DNAcopy.data, chr.col = "chrom", start.col = "loc.start", end.col = "loc.end", segment.value.col = "seg.mean", chr.transformation = "23:X,24:Y,25:MT", genome = NULL, verbose = TRUE)
#'
#' @param DNAcopy.data The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to "chrom")
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to "loc.start")
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to "loc.end")
#' @param segment.value.col The name or number of the column with segment value. If NULL, it is automatically identified. (default to "seg.mean")
#' @param chr.transformation  "23:X,24:Y,25:MT"
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#' 
#' library(DNAcopy)
#'
#' data(coriell)
#' CNA.object <- CNA(cbind(coriell$Coriell.05296), coriell$Chromosome, coriell$Position, data.type="logratio",sampleid="c05296")
#'
#' smoothed.CNA.object <- smooth.CNA(CNA.object)
#' segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
#'
#' DNAcopy.data <- segment.smoothed.CNA.object
#' loadCopyNumberCalls.DNAcopy(DNAcopy.data = DNAcopy.data)
#'
#'
#' @export loadCopyNumberCalls.DNAcopy
#'
#' @importFrom GenomicRanges mcols
#' 
loadCopyNumberCalls.DNAcopy <- function(DNAcopy.data,
                                        chr.col = "chrom", 
                                        start.col = "loc.start",
                                        end.col = "loc.end", 
                                        segment.value.col = "seg.mean", 
                                        chr.transformation = "23:X,24:Y,25:MT",
                                        genome = NULL,
                                        verbose = TRUE) {
  #if the Data come from DNAcopy
  if(methods::is(DNAcopy.data, "DNAcopy")) cnv.data <- data.frame(DNAcopy.data$output)
  
  cnv.data <- transformChr(seg.data = cnv.data, chr.transformation = chr.transformation)
  
  segs <- loadCopyNumberCalls(cnv.data = cnv.data, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = NA, 
                              segment.value.col = segment.value.col,
                              loh.col = NA, 
                              genome = genome,
                              verbose = verbose)
  

 

  return(segs)
}


