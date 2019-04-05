#' loadSNPData
#'
#' @description
#' Loads SNP array data in a tabular format
#'
#' @details
#' Given a file name or data in a tabular format,
#' the function loads SNP array data in a tabular format.
#' It will try to identify the columns with the relevant information
#' (chr, position, BAF, LRR, etc...) or will use the column number or name
#' supplied by the user, if any. It will convert the tabular data into a
#' GRanges, with one range per SNP in the table.
#'
#' @usage loadSNPData(snp.data, chr.col = NULL, start.col = NULL, end.col = NULL, pos.col = NULL, baf.col = NULL, lrr.col = NULL, id.col = NULL, genome = NULL, verbose = TRUE)
#'
#' @param snp.data Either the name of the file with the data or a variable containing the data.
#' @param start.col (number or character) the name or number of the column with start position
#' @param end.col (number or character) the name of number of the column with end position
#' @param genome (character) The name of the genome (defaults to "hg19")
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param pos.col (number or character) The name or number of the column with position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param baf.col (number or character) The name or number of the column with BAF information. If NULL, it is automatically identified. (defaults to NULL)
#' @param lrr.col (number or character) The name or number of the column with LRR information. If NULL, it is automatically identified. (defaults to NULL)
#' @param id.col (number or character) The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (defaults to NULL)
#' @param verbose Wether information messages should be generated. (defaults to TRUE)
#'
#' @return
#' A GRanges with a range per SNP
#'
#' @examples
#' 
#' # There are two examples of posible files to load.
#' snp.data1 <- system.file("extdata", "snp.data_test.csv", package = "CopyNumberPlots", mustWork=TRUE)
#' snps <- loadSNPData(snp.data = snp.data1)
#' 
#' snp.data2 <-  system.file("extdata", "snp.data_test2.csv", package = "CopyNumberPlots", mustWork=TRUE)
#' snps <- loadSNPData(snp.data = snp.data2)
#' 
#' @export loadSNPData
#' @importFrom utils read.table
#' @importFrom GenomicRanges mcols
#'
loadSNPData <- function(snp.data, 
                        chr.col = NULL,
                        start.col = NULL,
                        end.col = NULL, 
                        pos.col = NULL,
                        baf.col = NULL,
                        lrr.col = NULL, 
                        id.col = NULL,
                        genome = NULL,
                        verbose = TRUE) {
  
  #If its a file, try to load it
  if(is.character(snp.data)) {
    #if file dose not exist
    if(!file.exists(snp.data)){
      stop("The file ", snp.data, " does not exist")
    }
    
    if(verbose) message("Reading data from ", snp.data)
    
    #It is problematic with snparray data
    # #Try to load it using toGRanges
    # snps <- tryCatch(regioneR::toGRanges(snp.data, genome = genome), error = function(e){return(NULL)}, warning = function(w){})
    snps <- NULL
    
    #if the toGRanges failed, try with a series of read.tabe statements
    if(is.null(snps)) snps <- tryCatch(utils::read.table(snp.data, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    if(is.null(snps)) snps <- tryCatch(utils::read.table(snp.data, sep = ";", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    if(is.null(snps)) snps <- tryCatch(utils::read.table(snp.data, sep = ",", header = TRUE, stringsAsFactors = FALSE),
                                       error = function(e) return(NULL))
    if(is.null(snps)) stop("it was not possible to load the file")
    
    snp.data <- snps
  }
  
  #At this point, snp.data is either a GRanges or a data.frame
  
  #If it's  not a GRanges, try to convert it into a GRanges
  
  if(!methods::is(snp.data, "GRanges")) {
    #snps <- tryCatch(regioneR::toGRanges(snp.data, genome = genome), error = function(e){return(NULL)}, warning = function(w){})
    snps <- NULL
    if(is.null(snps)) { #If toGRanges failed try to identify columns by name
      chr.col <- getChrColumn(col = chr.col, df = snp.data, needed = TRUE,  verbose = verbose)
      start.col <- getStartColumn(col = start.col, df = snp.data, needed = FALSE, verbose = verbose)
      end.col <- getEndColumn(col = end.col, df = snp.data, needed = FALSE, verbose = verbose)
      pos.col <- getPosColumn(col = pos.col, df = snp.data, needed = FALSE, verbose = verbose)
      
     
      # chr.col
      if(!is.null(chr.col)) names(snp.data)[chr.col] <- "chr"
      
      # snp position
      
      if (!is.null(pos.col)) {
        start.col <- pos.col
        end.col <- pos.col
      } else {
        if(is.null(start.col) || is.null(end.col)){
          stop("It was not possible to identify the required data: either Start and End or Position is requierd")
        } 
      }  
      

      other.cols <- seq_len(length(snp.data))[!(seq_len(length(snp.data)) %in% c(chr.col, start.col, end.col))]
      columns <- c(chr.col, start.col, end.col, other.cols)
      snps <- tryCatch(regioneR::toGRanges(snp.data[,columns], genome = genome), 
                       error = function(e){
                         stop("It was not possible to transform the data into a GRanges. Is there any format specific data loading function available? ", e)
                       },
                       warning = function(w){})
    }
    if(!methods::is(snps, "GRanges")) stop("It was not possible to read and transform the data. Is there any format specific data loading function available?")
  } else {
    #If it is a GRanges, simply change its name
    snps <- snp.data
  }
  
  #BAF
  baf.col <- getBAFColumn(df = GenomicRanges::mcols(snps), col = baf.col, needed = FALSE, verbose = verbose)
  if(!is.null(baf.col))  names(GenomicRanges::mcols(snps))[baf.col] <- "baf"

  
  #LRR
  lrr.col <- getLRRColumn(df = GenomicRanges::mcols(snps), col = lrr.col, needed = FALSE, verbose = verbose)
  if(!is.null(lrr.col)) names(GenomicRanges::mcols(snps))[lrr.col] <- "lrr"
  
  #ID
  id.col <-  getIDColumn(df = GenomicRanges::mcols(snps), col = id.col, avoid.pattern = "sample", needed = FALSE, verbose = verbose)
  if(!is.null(id.col))names(GenomicRanges::mcols(snps))[id.col] <- "id"
  
  snps <- sort(snps)
  
  return(snps)
}


