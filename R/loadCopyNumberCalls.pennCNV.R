#' loadCopyNumberCalls.pennCNV
#'
#' @description
#' Loads copy number calls from pennCNV.rawcnv file format 
#'
#' @details
#' This function will load copy number calls from pennCNV.rawcnv file format. 
#' Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files. 
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#' @usage loadCopyNumberCalls.pennCNV (pennCNV.file, chr.col = NULL, start.col = NULL, end.col = NULL, cn.col = NULL, segment.value.col = NULL, genome = NULL, verbose=TRUE)
#'
#' @param pennCNV.file The name of the file with the data
#' @param chr.col (number or character) The name or number of the column with chromosome information. If NULL, it is automatically identified. (defaults to NULL)
#' @param start.col (number or character) The name or number of the column with start position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param end.col (number or character) The name or number of the column with end position information. If NULL, it is automatically identified. (defaults to NULL)
#' @param segment.value.col (number or character) The name or number of the column with segment value. If NULL, it is automatically identified. (defaults to NULL)
#' @param cn.col (number or character) The name or number of the column with CN information. If NULL, it is automatically identified. (defaults to NULL)
#' @param genome (character) The name of the genome (defaults to NULL)
#' @param verbose (logical) Whether to show information messages. (defaults to TRUE)
#' 
#' @return
#' A GRanges with a range per copy number segment or a list of GRanges with a GRanges per sample.
#'
#' @examples
#' ## loadCopyNumberCalls.pennCNV from .rawcnv file format: 
#' ## The file to run the example can be found in: http://penncnv.openbioinformatics.org/en/latest/user-guide/test/
#' 
#' pennCNV.file <- system.file("extdata", "pennCNV.rawcnv", package = "CopyNumberPlots", mustWork = TRUE)
#' cnv.call <- loadCopyNumberCalls.pennCNV(pennCNV.file = pennCNV.file)
#'
#' @export loadCopyNumberCalls.pennCNV
#'

loadCopyNumberCalls.pennCNV <- function(pennCNV.file,
                                        chr.col = NULL, 
                                        start.col = NULL,
                                        end.col = NULL, 
                                        cn.col = NULL,
                                        segment.value.col = NULL, 
                                        genome = NULL, 
                                        verbose=TRUE){
  
  penncnv.info <- utils::read.table(file = pennCNV.file, stringsAsFactors = FALSE)
  
  # variables for chr star and end
  penncnv.data <- regioneR::toDataframe(regioneR::toGRanges(penncnv.info$V1))
  
  #variables for numsnp
  numsnp <- do.call(rbind, strsplit(x = penncnv.info$V2, split = "="))
  
  if(length(unique(numsnp[,1])) != 1) stop("There are more than one variable in numsnp column")
  
  df <- data.frame(as.numeric(numsnp[,2]))
  names(df) <- unique(numsnp[,1])
  penncnv.data <- cbind(penncnv.data, df)
  
  #variable for lengthsnp
  lengthsnp <- do.call(rbind, strsplit(x = penncnv.info$V3, split = "="))
  lengthsnp[,2] <- gsub(pattern = ",", replacement = "", x = lengthsnp[,2])
  df <- data.frame(lengthsnp[,2])
  
  if(length(unique(lengthsnp[,1])) != 1) stop("There are more than one variable in length column")
  
  df <- data.frame(as.numeric(lengthsnp[,2]))
  names(df) <- unique(lengthsnp[,1])
  penncnv.data <- cbind(penncnv.data, df)
  
  
  #variable for states and cn
  states.cn <- do.call(rbind, strsplit(x=penncnv.info$V4, split = ","))
  states <- states.cn[,1] 
  penncnv.data <- cbind(penncnv.data, states)
  cn <- do.call(rbind, strsplit(x= states.cn[,2], split = "="))
  
  if(length(unique(cn[,1])) != 1) stop("There are more than one variable in cn column")
  
  df <- data.frame(as.numeric(cn[,2]))
  names(df) <- unique(cn[,1])
  penncnv.data <- cbind(penncnv.data, df)
  
  
  #variable for samples
  samples <- gsub(pattern = ".txt", replacement = "", x = penncnv.info$V5)
  penncnv.data <- cbind(penncnv.data, samples)
  
  #variable for startsnp
  startsnp <-  do.call(rbind, strsplit(x = penncnv.info$V6, split = "="))
  
  if(length(unique(startsnp[,1])) != 1) stop("There are more than one variable in startsnp column")
  
  df <- data.frame(startsnp[,2])
  names(df) <- unique(startsnp[,1])
  penncnv.data <- cbind(penncnv.data, df)
  
  
  #variable for endsnp
  endsnp <-  do.call(rbind, strsplit(x = penncnv.info$V7, split = "="))
  
  if(length(unique(endsnp[,1])) != 1) stop("There are more than one variable in endsnp column")
  
  df <- data.frame(endsnp[,2])
  names(df) <- unique(endsnp[,1])
  penncnv.data <- cbind(penncnv.data, df)
  
  
  segs <- loadCopyNumberCalls(cnv.data = penncnv.data, 
                              chr.col = "chr", 
                              start.col = "start",
                              end.col = "end", 
                              cn.col = "cn",
                              segment.value.col = NULL, 
                              genome = genome, 
                              verbose=TRUE)
  
  
  #if there is more than one sample in pennCNV results
  if(length(table(segs$samples)) > 1){
    segs <- split(segs, GenomicRanges::mcols(segs)$samples)
  }
  
  return(segs)
}