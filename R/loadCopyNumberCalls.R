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
#' @usage loadCopyNumberCalls(cnv.data, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = NULL, loh.col = NULL, segment.value.col = NULL, genome = NULL, zero.based = FALSE, verbose = TRUE)
#'
#' @param cnv.data Either the name of the file with the data or a variable containing the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param loh.col (number or character) The name or number of the column with LOH information. If NULL, it is automatically identified. (defaults to NULL)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to NULL)
#' @param zero.based (logical) Whether the data is zero-based and half open (i.e. ranges are defined by (start:end] so chr1:10-20 represents nine bases long features spanning from base 11 to 20). (defaults to FALSE)
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per copy number segment
#'
#' @examples
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' 
#' cnv.call <- loadCopyNumberCalls(cnv.data = df) 
#'
#' @export loadCopyNumberCalls
#'
#' @importFrom GenomicRanges mcols start
#' 

loadCopyNumberCalls <- function(cnv.data,
                                chr.col = NULL,
                                start.col = NULL, 
                                end.col = NULL, 
                                cn.col = NULL,
                                loh.col = NULL,
                                segment.value.col = NULL,
                                genome = NULL,
                                zero.based = FALSE,
                                verbose = TRUE){
  
  #If its a file, try to load it
  if(is.character(cnv.data)) {
    #if file dose not exist
    if(!file.exists(cnv.data)){
      stop("The file ", cnv.data, " does not exist")
    }
    
    if(verbose) message("Reading data from ", cnv.data)
    
    #if file exists try to load it using toGRanges
    segs <- tryCatch(regioneR::toGRanges(cnv.data), error = function(e){return(NULL)}, warning = function(w){})
    #if the toGRanges failed, try with a series of read.tabe statements
    if(is.null(segs)) segs <- tryCatch(utils::read.table(cnv.data, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    if(is.null(segs)) segs <- tryCatch(utils::read.table(cnv.data, sep = ";", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    if(is.null(segs)) segs <- tryCatch(utils::read.table(cnv.data, sep = ",", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    
    if(is.null(segs)) stop("it was not possible to load the file")
    
    cnv.data <- segs
  }
  
  #At this point, cnv.data is either a GRanges or a data.frame
  
  #If it's  not a GRanges, try to convert it into a GRanges
  if(!methods::is(cnv.data, "GRanges")) {
    segs <- tryCatch(regioneR::toGRanges(cnv.data, genome = genome), error = function(e){}, warning = function(w){})
    if(is.null(segs)) { #If toGRanges failed
      chr.col <- getChrColumn(col = chr.col, df = cnv.data, needed = TRUE, verbose = verbose)
      start.col <- getStartColumn(col = start.col, df = cnv.data, needed = TRUE, verbose = verbose)
      end.col <- getEndColumn(col = end.col, df = cnv.data, needed = TRUE, verbose = verbose)
      
      #This if should never be true. It should have failed in the column identification functions
      if(any(is.null(chr.col), is.null(start.col), is.null(end.col))) stop("It was not possible to identify the required data: Chromosome, Start and End")
      
      other.cols <- seq_len(length(cnv.data))[-c(chr.col, start.col, end.col)]
      segs <- tryCatch(regioneR::toGRanges(cnv.data[,c(chr.col, start.col, end.col, other.cols)], genome = genome), 
                       error = function(e){
                         stop("It was not possible to transform the data into a GRanges. Is there any format specific data loading function available? ", e)
                       },
                       warning = function(w){})
    }
    if(!methods::is(segs, "GRanges")) stop("It was not possible to read and transform the data. Is there any format specific data loading function available?")
  } else {
    #If it is a GRanges, simply change its name
    segs <- cnv.data
  }
  
  #If we are here, we have a GRanges with our data
    #identify the columns we want
    #Copy Number
    cn.col <- getCopyNumberColumn(df = GenomicRanges::mcols(segs), col = cn.col, needed = FALSE, verbose = verbose)
    if(!is.null(cn.col)) names(GenomicRanges::mcols(segs))[cn.col] <- "cn"

    #LOH
    loh.col <- getLOHColumn(df = GenomicRanges::mcols(segs), col = loh.col, needed = FALSE, verbose = verbose)
    if(!is.null(loh.col)) names(GenomicRanges::mcols(segs))[loh.col] <- "loh"

    #Segment Value
    segment.value.col <- getSegmentValueColumn(df = GenomicRanges::mcols(segs), col = segment.value.col, needed = FALSE, verbose = verbose)
    if(!is.null(segment.value.col)) names(GenomicRanges::mcols(segs))[segment.value.col] <- "segment.value"
    
    #If segs ranges of data is zero.based
    if(isTRUE(zero.based)) GenomicRanges::start(segs)<- GenomicRanges::start(segs)+1
    

  return(segs)
}
