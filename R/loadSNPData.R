#' loadSNPData
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
#' @usage loadSNPData(snps.file, genome="hg19", chr.col=NULL, pos.col=NULL, baf.col=NULL, lrr.col=NULL, snp.col=NULL, verbose=TRUE)
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
#' @export loadSNPData
#' @importFrom regioneR toGRanges
#' @importFrom utils read.table
#'

#Read the SNPs raw data file, identify the relevant columns
loadSNPData <- function(snps.file, genome="hg19", chr.col=NULL, pos.col=NULL, baf.col=NULL, lrr.col=NULL, snp.col=NULL, verbose=TRUE) {
  if(verbose==TRUE) message("Reading file ", snps.file, "...")
  snps <- read.table(snps.file, sep="\t", header=TRUE, stringsAsFactors = FALSE)

  #Get the columns
  if(verbose==TRUE) message("Indentifying data columns...")
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
        stop("The column ", msg.col.name, " was not found in the data file")
      } else {
        col.num <- integer(0)
      }
    }
    return(col.num)
  }

  chr.col <- getColumn(chr.col, "Chr|chr", snps, "Chromosome" )
  pos.col <- getColumn(pos.col, "Position|Pos|pos", snps, "Position" )
  baf.col <- getColumn(baf.col, "BAF|B.Allele|Freq|freq", snps, "B-Allele Frequency" )
  lrr.col <- getColumn(lrr.col, "LRR|Log.R.Ratio|Log", snps, "Log Ratio")
  snp.col <- getColumn(snp.col, "Name|Id|SNP.Name|SNP.id", snps, "SNP Identifier", needed=FALSE)

  if(length(snp.col)>0) {
    snps <- snps[,c(chr.col, pos.col, pos.col, baf.col, lrr.col, snp.col)]
    names(snps) <- c("chr", "start", "end", "baf", "lrr", "snp.name")
  } else {
    snps <- snps[,c(chr.col, pos.col, pos.col, baf.col, lrr.col)]
    names(snps) <- c("chr", "start", "end", "baf", "lrr")
  }
  snps <- regioneR::toGRanges(snps)
  seqlevelsStyle(snps) <- "UCSC"

  if(verbose==TRUE) message("Removing SNPs with NA values...")
  baf.is.na <- is.na(snps$baf)
  lrr.is.na <- is.na(snps$lrr)
  if(verbose==TRUE) message("BAF is NA: ", length(which(baf.is.na)), "   LRR is NA: ", length(which(lrr.is.na)), "   Total: ", length(which(baf.is.na | lrr.is.na)))
  snps <- snps[!(baf.is.na | lrr.is.na)]

  if(verbose==TRUE) message("Removing SNPs out of the canonical chromosomes...")
  length.before <- length(snps)
  snps <- filterChromosomes(snps, organism = genome, chr.type = "canonical")
  if(verbose==TRUE) message("Removed ", length.before - length(snps), " SNPs in non-canonical chromosomes")

  snps <- sort(snps)

  return(snps)
}
