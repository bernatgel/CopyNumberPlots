#' loadSNPDataFromVCF
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
#' @usage loadSNPDataFromVCF(vcf.file, genome="hg19", randomize.baf=TRUE, verbose=TRUE)
#'
#' @param snps.file The name of the file with the data
#' @param genome The name of the genome (default to "hg19")
#' @param randomize.baf Flip the baf of about half the snps (the ones in odd positions in the genome) to achieve a mirror-like effect as in SNP arrays (defaults to TRUE)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per SNP
#'
#' @examples
#'
#'
#'
#' @export loadSNPDataFromVCF
#' @import VariantAnnotation

#Read a VCF file and extract coverage (LRR) and frequency (BAF) information from it
loadSNPDataFromVCF <- function(vcf.file, genome="hg19", randomize.baf=TRUE, verbose=TRUE) {
  if(verbose==TRUE) message("Scanning file ", vcf.file, "...")

  vcf.header <- scanVcfHeader(vcf.file)
  if(!("AD" %in% row.names(geno(vcf.header)))) stop("The VCF file does not have the AD field in genotype. BAF/LRR computation from FREQ and DP still not implemented.")

  #TODO: Should we accept a GRanges to scan only specific regions of the file? (for zooming, etc...)
  vars <- readVcf(file=TabixFile(vcf.file), genome = "hg19", param = ScanVcfParam(info=NA, geno = "AD"))

  #TODO: Remove the indels?

  #BAF
    #Compute the freq
    ad <- geno(vars)$AD
    ad.ref <- unlist(lapply(ad, "[", 1))
    ad.alt <- unlist(lapply(ad, "[", 2))

    baf <- ad.alt/(ad.alt+ad.ref)

  #LRR
    lrr <- ad.alt+ad.ref+1
    lrr <- log(lrr/mean(lrr))

  #Build the GRanges
  vars <- rowRanges(vars)
  mcols(vars) <- data.frame(baf=baf, lrr=lrr)

  if(randomize.baf==TRUE) {
    #flip the frequency of half the SNPs to achieve the same mirroring option we see in SNP-arrays
    odd.pos <- start(vars) %% 2 == 1
    vars[odd.pos]$baf <- 1 - vars[odd.pos]$baf
  }

  vars <- sort(vars)

  return(vars)
}


