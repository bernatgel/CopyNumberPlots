#' loadCopyNumberCallsDNAcopy
#'
#' @description
#' Loads copy number calls from DNACopy results.
#'
#' @details
#' This function will load segments data from DNAcopy data structure
#' in R. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files.
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCallsDNAcopy(DNAcopy.data, chr.col = "chrom", start.col = "loc.start", end.col = "loc.end", segment.value.col = "seg.mean", cn.col = NULL, chr.transformation = "23:X,24:Y,25:MT", genome = NULL, verbose = TRUE)
#'
#' @param DNAcopy.data The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to "chrom")
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to "loc.start")
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to "loc.end")
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to "seg.mean")
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param chr.transformation (character) The transformation of the chromosome names in a comma separated "key:value" format.(defaults to "23:X,24:Y,25:MT")
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment or a list of GRanges with a GRanges per sample.
#'
#' @examples
#'
#' library(DNAcopy)
#'
#' data(coriell)
#' CNA.object <- suppressWarnings(CNA(cbind(coriell$Coriell.05296), coriell$Chromosome, coriell$Position, data.type="logratio",sampleid="c05296"))
#'
#' smoothed.CNA.object <- smooth.CNA(CNA.object)
#' DNAcopy.data <- segment(smoothed.CNA.object, verbose=1)
#'
#' cnv.call <- loadCopyNumberCallsDNAcopy(DNAcopy.data = DNAcopy.data)
#'
#' # more than 1 sample
#' CNA.object <- CNA(genomdat = cbind(coriell$Coriell.05296, coriell$Coriell.13330), chrom = coriell$Chromosome, maploc = coriell$Position, data.type = "logratio", sampleid = c("c05296", "c13330"))
#' smoothed.CNA.object <- smooth.CNA(CNA.object)
#' DNAcopy.data <- segment(smoothed.CNA.object, verbose=1)
#'
#' cnv.call <- loadCopyNumberCallsDNAcopy(DNAcopy.data = DNAcopy.data)
#'
#' @export loadCopyNumberCallsDNAcopy
#'
loadCopyNumberCallsDNAcopy <- function(DNAcopy.data,
                                        chr.col = "chrom",
                                        start.col = "loc.start",
                                        end.col = "loc.end",
                                        segment.value.col = "seg.mean",
                                        cn.col = NULL,
                                        chr.transformation = "23:X,24:Y,25:MT",
                                        genome = NULL,
                                        verbose = TRUE) {
  #if the Data come from DNAcopy
  if(methods::is(DNAcopy.data, "DNAcopy")){
    cnv.data <- data.frame(DNAcopy.data$output)
  }else{
    stop("DNAcopy.data must be DNAcopy object")
  }


  if(!is.null(chr.transformation)){
    chr.col <- names(cnv.data)[getChrColumn(df = cnv.data, col = chr.col, verbose = FALSE)]
    cnv.data[,chr.col] <- transformChr(chr = cnv.data[,chr.col], chr.transformation = chr.transformation)
  }


  segs <- loadCopyNumberCalls(cnv.data = cnv.data,
                              chr.col = chr.col,
                              start.col = start.col,
                              end.col = end.col,
                              cn.col = cn.col,
                              segment.value.col = segment.value.col,
                              loh.col = NULL,
                              genome = genome,
                              verbose = verbose)


  #if there is more than one sample in segs

  if(length(unique(segs$ID)) > 1){
    segs <- split(segs, mcols(segs)$ID)
  }


  return(segs)
}


