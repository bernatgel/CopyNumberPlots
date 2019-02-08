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
#' @usage loadCopyNumberCalls(cnv.data, chr.col=NULL, pos.col=NULL, cn.col=NULL, loh.col=NULL, segment.value.col=NULL, genome=NULL, verbose=TRUE)
#'
#' @param snps.file The name of the file with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param pos.col The name or number of the column with position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with BAF information. If NULL, it is automatically identified. (default to NULL)
#' @param loh.col The name or number of the column with LRR information. If NULL, it is automatically identified. (default to NULL)
#' @param segment.value.col The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (default to NULL)
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per SNP
#'
#' @examples
#'
#'
#'
#' @export loadCopyNumberCalls


#Function used in loadSNPData. Make it an internal function shared by both functions

getColumn <- function(col, pattern, df, msg.col.name, needed=TRUE) {
  col.num <- integer(0)
  if(is.null(col)) {
    col.num <- which(grepl(names(df), pattern = pattern))[1]
  } else {
    if(is.numeric(col)) {
      col.num <- col
    } else {
      if(is.character(col))  {
        col.num <- which(names(df)==col)
      }
    }
  }
  if(is.na(col.num) || length(col.num)==0) {
    if(needed==TRUE) {
      stop("The column ", msg.col.name, " was not found in the data")
    } else {
      col.num <- NULL
    }
  }
  return(col.num)
}

#return the column representing the chromosomes
getChrColumn <- function(chr.col, df, needed=TRUE) {
  return(getColumn(chr.col, "Chr|chr", df, "Chromosome", needed=needed))
}

getPosColumn <- function(pos.col, df, needed=TRUE) {
  return(getColumn(pos.col, "Position|Pos|pos|loc|maploc", df, "Position",needed=needed))
}

getStartColumn <- function(start.col, df, needed=TRUE) {
  return(getColumn(pos.col, "Start|start|First|first|Begin|begin", df, "Start",needed=needed))
}

getEndColumn <- function(end.col, df, needed=TRUE) {
  return(getColumn(end.col, "End|end|Last|last", df, "End", needed=needed))
}

getCopyNumberColumn <- function(cn.col, df, needed=TRUE) {
  return(getColumn(cn.col, "CN|cn|Copy|copy", df, "Copy Number", needed=needed))
}

getLOHColumn <- function(loh.col, df, needed=TRUE) {
  return(getColumn(cn.col, "LOH|Loh|loh|Loss|loss", df, "LOH", needed=needed))
}

getSegmentValueColumn <- function(segment.value.col, df, needed=TRUE) {
  return(getColumn(cn.col, "Value|value|mean|median|ratio", df, "Segment Value", needed=needed))
}




loadCopyNumberCalls <- function(cnv.data, chr.col=NULL, start.col=NULL, end.col=NULL, cn.col=NULL, loh.col=NULL, segment.value.col=NULL, genome=NULL, verbose=TRUE) {
  #If its a file, try to load it
  if(is.character(cnv.data)) {
    if(verbose) message("Reading data from ", cnv.data)
    #try to load it using toGRanges
    #if failed
      #try to load it using read.table
      #Should we accept a sep parameter? and a quote? skip? No, read the table and process the data.frame
      #if failed => error
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
      cn.col <- getCopyNumberColumn(cn.col, mcols(segs), needed=FALSE)
      if(!is.null(cn.col)) names(mcols(segs))[cn.col] <- "cn"
    #LOH
      loh.col <- getLOHColumn(loh.col, mcols(segs), needed=FALSE)
      if(!is.null(loh.col)) names(mcols(segs))[loh.col] <- "loh"
    #Segment Value
      segment.value.col <- getSegmentValueColumn(segment.value.col, mcols(segs), needed = FALSE)
      if(!is.null(segment.value.col)) names(mcols(segs))[segment.value.col] <- "segment.value"

  return(segs)
}

Fer una funcio per a fer 2n el que no tingui segments. Afegir com a pas extra aqui? Hauria de deixar totes les columnes no conegudes com a NA

Les funcions específiques de paquet podrien dependre d aquesta amb paràmetres prefixats per a reduïr codi i augmentar funcionalitat

loadCopyNumberCalls(cnv.data)
