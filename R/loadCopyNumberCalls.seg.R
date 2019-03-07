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
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls.seg(seg.file, chr.col = "chrom", start.col = "loc.start", end.col = "loc.end", segment.value.col = "seg.mean", chr.transformation = "23:X,24:Y,25:MT", genome = NULL, verbose = TRUE)
#'
#' @param seg.file The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to "chrom")
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to "loc.start")
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to "loc.end")
#' @param segment.value.col The name or number of the column with segment value. If NULL, it is automatically identified. (default to "seg.mean")
#' @param chr.transformation Transformation of the chromosomes. The files.seg have the chromosomes in a numerical form.
#' In order to transform the sexual and mitocondrial chromosomes from numbers to characters,
#' the user has to provide in the chr.transformation parameter as describe in 
#' https://cnvkit.readthedocs.io/en/stable/importexport.html#import-seg, in a python dictionary format ("key:value").(default to "23:X,24:Y,25:MT" )
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#' ## loadCopyNumberCalls.seg from .seg file format: the file to run the example can be found in: https://software.broadinstitute.org/software/igv/SEG
#' ## under example.seg file name.
#' seg.file <- "testFiles/example.seg"
#' loadCopyNumberCalls.seg(seg.file)
#'
#' @export loadCopyNumberCalls.seg
#'

loadCopyNumberCalls.seg <- function(seg.file, 
                                    chr.col = "chrom", 
                                    start.col = "loc.start",
                                    end.col = "loc.end", 
                                    segment.value.col = "seg.mean", 
                                    genome = NULL,
                                    chr.transformation = "23:X,24:Y,25:MT",
                                    verbose = TRUE){
  
  #Load file.seg. First, probe if file exist
  if(!file.exists(seg.file)) stop(paste0(seg.file, " does not exist or you are not in the correct directory."))
  
  seg.data <- read.table(file = seg.file, sep = "\t", skip = 1, stringsAsFactors = FALSE)
  colnames(seg.data) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  
  # we change the order of the columns to be transfromed in a GRanges
  seg.data <- seg.data[,c("chrom", "loc.start", "loc.end", "ID", "num.mark", "seg.mean")] 
    
  # Transformation of the chromosomes. The segs files have the chromosomes in a numerical form
  # In order to transform the sexual and mitocondrial chromosomes numbers in characters,
  # we transform them using this funciton and giving a chr.transformation parameter in the same
  # way than in the example to transform it.
  seg.data <- transformChr(seg.data = seg.data, chr.transformation = chr.transformation)
    
  segs <- loadCopyNumberCalls(cnv.data = seg.data, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = NA, 
                              segment.value.col = segment.value.col,
                              loh.col = NA, 
                              genome = genome,
                              verbose = verbose)
  
  # Once we have a GRanges, transform the chromosomes in the UCSC style
  GenomeInfoDb::seqlevelsStyle(segs) <- "UCSC"
  
  #if there is more than one sample in cn.mops.res
  if(length(table(segs$ID)) > 1){
    segs <- split(segs, mcols(segs)$ID)
  }
  
  return(segs)
}
