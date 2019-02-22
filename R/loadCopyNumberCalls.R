#' loadCopyNumberCalls
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
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls(cnv.data, chr.col=NULL, start.col=NULL, end.col=NULL, cn.col=NULL, loh.col=NULL, segment.value.col=NULL, genome=NULL, verbose=TRUE)
#'
#' @param cnv.data The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to NULL)
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with BAF information. If NULL, it is automatically identified. (default to NULL)
#' @param loh.col The name or number of the column with LRR information. If NULL, it is automatically identified. (default to NULL)
#' @param segment.value.col The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (default to NULL)
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#'
#'
#'
#' @export loadCopyNumberCalls
#'
#' @importFrom GenomicRanges mcols


#IDEA: create a function to make every part of the genome not covered by  the segments, a 2n segment. Do NOT call it from here automagically.

loadCopyNumberCalls <- function(cnv.data, chr.col=NULL, start.col=NULL, end.col=NULL, cn.col=NULL, loh.col=NULL, segment.value.col=NULL, genome=NULL, verbose=TRUE) {
  #If its a file, try to load it
  if(is.character(cnv.data)) {
    if(verbose) message("Reading data from ", cnv.data)
    
    #if file dose not exist
    if(!file.exists(cnv.data)){
      stop(paste0(cnv.data, " not exist"))
    }
    
    #if file exists try to load it using toGRanges
    gr <- tryCatch(toGRanges(cnv.data, genome=genome),error = function(e){return(NULL)})
    #if failed
    if(is.null(gr)){
      #try to load it using read.table
      #Should we accept a sep parameter? and a quote? skip? No, read the table and process the data.frame
      #if failed => error
      stop("cannot load the file. Please give a GRanges object class  to read cnv.data")
    } else{
      segs <- tryCatch(toGRanges(cnv.data, genome=genome), error=function(e) {})
      if(is.null(segs)) { #If toGRanges failed
        chr.col <- getChrColumn(chr.col, cnv.data, needed=TRUE)
        start.col <- getStartColumn(start.col, cnv.data, needed=TRUE)
        end.col <- getEndColumn(end.col, cnv.data, needed=TRUE)

        #This if should never be true. It should have failed in the column identification functions
        if(any(is.null(chr.col), is.null(start.col), is.null(end.col))) stop("It was not possible to identify the required data: Chromosome, Start and End")

        other.cols <- seq_len(length(cnv.data))[-c(chr.col, start.col, end.col)][-start.col][-end.col]
        segs <- toGRanges(cnv.data[,c(chr.col, start.col, end.col, other.cols)], genome=genome)
      }
      if(!methods::is(segs, "GRanges")) stop("It was not possible to read and transform the data. Is there any format specific data loading function available?")

    }
    
  }

  #If it's  not a GRanges, try to convert it into a GRanges
  if(!methods::is(cnv.data, "GRanges")) {
    segs <- tryCatch(toGRanges(cnv.data, genome=genome), error=function(e) {})
    if(is.null(segs)) { #If toGRanges failed
      chr.col <- getChrColumn(chr.col, cnv.data, needed=TRUE)
      start.col <- getStartColumn(start.col, cnv.data, needed=TRUE)
      end.col <- getEndColumn(end.col, cnv.data, needed=TRUE)

      #This if should never be true. It should have failed in the column identification functions
      if(any(is.null(chr.col), is.null(start.col), is.null(end.col))) stop("It was not possible to identify the required data: Chromosome, Start and End")

      other.cols <- seq_len(length(cnv.data))[-c(chr.col, start.col, end.col)][-start.col][-end.col]
      segs <- toGRanges(cnv.data[,c(chr.col, start.col, end.col, other.cols)], genome=genome)
    }
    if(!methods::is(segs, "GRanges")) stop("It was not possible to read and transform the data. Is there any format specific data loading function available?")
  }

  #If we are here, we have a GRanges with our data
    #identify the columns we want
    #Copy Number
      cn.col <- getCopyNumberColumn(GenomicRanges::mcols(segs), cn.col,needed=FALSE)
      if(!is.null(cn.col)) names(GenomicRanges::mcols(segs))[cn.col] <- "cn"
    #LOH
      loh.col <- getLOHColumn(GenomicRanges::mcols(segs), loh.col, needed=FALSE)
      if(!is.null(loh.col)) names(GenomicRanges::mcols(segs))[loh.col] <- "loh"
    #Segment Value
      segment.value.col <- getSegmentValueColumn(GenomicRanges::mcols(segs), segment.value.col, needed = FALSE)
      if(!is.null(segment.value.col)) names(GenomicRanges::mcols(segs))[segment.value.col] <- "segment.value"

  return(segs)
}
