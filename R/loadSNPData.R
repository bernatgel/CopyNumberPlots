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
#'snp.file <- "test.snp.file.csv"
#'
#'snp.df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#'                     "copy.number.level"=3, "LOH"=0,"median.value.per.segment"=1.2, 
#'                     "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#'
#'write.table(x = snp.df, file = snp.file, col.names = TRUE, row.names = FALSE, sep = "\t")
#'
#'loadSNPData(snps.file = snp.file)
#'
#'
#' @export loadSNPData
#' @importFrom regioneR toGRanges filterChromosomes
#' @importFrom utils read.table
#' @importFrom GenomeInfoDb seqlevelsStyle
#'

#Read the SNPs raw data file, identify the relevant columns
loadSNPData <- function(snps.file, genome = "hg19", chr.col = NULL, pos.col = NULL, baf.col = NULL, lrr.col = NULL, id.col = NULL, verbose = TRUE) {
  if(verbose==TRUE) message("Reading file ", snps.file, "...")
  snps <- read.table(snps.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  #Get the columns
  chr.col <- getChrColumn(df = snps, col = chr.col, needed = TRUE)
  pos.col <- getPosColumn(df = snps, col = pos.col, needed = TRUE)
  baf.col <- getBAFColumn(df = snps, col = baf.col, needed = TRUE)
  lrr.col <- getLRRColumn(df = snps, col = lrr.col, needed = TRUE)
  id.col <- getIDColumn(df = snps, col = id.col, needed = TRUE)
  
 
  if(length(id.col)>0) {
    snps <- snps[,c(chr.col, pos.col, pos.col, baf.col, lrr.col, id.col)]
    names(snps) <- c("chr", "start", "end", "baf", "lrr", "id.name")
  } else {
    snps <- snps[,c(chr.col, pos.col, pos.col, baf.col, lrr.col)]
    names(snps) <- c("chr", "start", "end", "baf", "lrr")
  }
  snps <- regioneR::toGRanges(snps)
  GenomeInfoDb::seqlevelsStyle(snps) <- "UCSC"

  if(verbose==TRUE) message("Removing SNPs with NA values...")
  baf.is.na <- is.na(snps$baf)
  lrr.is.na <- is.na(snps$lrr)
  if(verbose==TRUE) message("BAF is NA: ", length(which(baf.is.na)), "   LRR is NA: ", length(which(lrr.is.na)), "   Total: ", length(which(baf.is.na | lrr.is.na)))
  snps <- snps[!(baf.is.na | lrr.is.na)]

  if(verbose==TRUE) message("Removing SNPs out of the canonical chromosomes...")
  length.before <- length(snps)
  snps <- regioneR::filterChromosomes(snps, organism = genome, chr.type = "canonical")
  if(verbose==TRUE) message("Removed ", length.before - length(snps), " SNPs in non-canonical chromosomes")

  snps <- sort(snps)

  #TODO: Set it to the original seqlevelsStyle?

  return(snps)
}
