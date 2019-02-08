#' loadCopyNumberCalls
#'
#' @description
#' Loads SNP array data in a tabular format
#'
#' @details
#' Given a file name, the function loads SNP array data in a tabular format.
#' It will try to identify the columns with the relevant information
#' (chr, position, BAF, LRR, etc...) or will use the column number or name
#' supplied by the user, if any. It will convert the tabular data into a
#' GRanges, with one range per SNP in the table.
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls(snps.file, genome="hg19", chr.col=NULL, pos.col=NULL, baf.col=NULL, lrr.col=NULL, snp.col=NULL, verbose=TRUE)
#'
#' @param snps.file The name of the file with the data
#' @param genome The name of the genome (default to "hg19")
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param pos.col The name or number of the column with position information. If NULL, it is automatically identified. (default to NULL)
#' @param baf.col The name or number of the column with BAF information. If NULL, it is automatically identified. (default to NULL)
#' @param lrr.col The name or number of the column with LRR information. If NULL, it is automatically identified. (default to NULL)
#' @param snp.col The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (default to NULL)
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



loadCopyNumberCalls <- function(snps.file, genome="hg19", chr.col=NULL, pos.col=NULL, baf.col=NULL, lrr.col=NULL, snp.col=NULL, verbose=TRUE) {

}
