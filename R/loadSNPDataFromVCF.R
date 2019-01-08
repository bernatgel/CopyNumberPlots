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
#' @param vcf.file The name of the file with the data
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
#' @importFrom VariantAnnotation readVcf ScanVcfParam scanVcfHeader geno
#' @importFrom Rsamtools TabixFile
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges start

#Read a VCF file and extract coverage (LRR) and frequency (BAF) information from it
loadSNPDataFromVCF <- function(vcf.file, genome="hg19", randomize.baf=TRUE, verbose=TRUE) {
  if(verbose==TRUE) message("Scanning file ", vcf.file, "...")

  vcf.header <- VariantAnnotation::scanVcfHeader(vcf.file)
  if(!("AD" %in% row.names(VariantAnnotation::geno(vcf.header)))) stop("The VCF file does not have the AD field in genotype. BAF/LRR computation from FREQ and DP still not implemented.")

  #TODO: Should we accept a GRanges to scan only specific regions of the file? (for zooming, etc...)
  vars <- VariantAnnotation::readVcf(file=Rsamtools::TabixFile(vcf.file), genome = "hg19", param = VariantAnnotation::ScanVcfParam(info=NA, geno = "AD"))

  #TODO: Remove the indels? Or at least make sure the frequency is correct

  #BAF
    #Compute the freq
    ad <- VariantAnnotation::geno(vars)$AD
    ad.ref <- unlist(lapply(ad, "[", 1))
    ad.alt <- unlist(lapply(ad, "[", 2))

    baf <- ad.alt/(ad.alt+ad.ref)

  #LRR
    lrr <- ad.alt+ad.ref+1
    lrr <- log(lrr/mean(lrr))

  #Build the GRanges
  vars <- SummarizedExperiment::rowRanges(vars)
  GenomicRanges::mcols(vars) <- data.frame(baf=baf, lrr=lrr)

  if(randomize.baf==TRUE) {
    #flip the frequency of every other SNP to achieve the same mirroring effect
    #we see in SNP-arrays
    odd.pos <- GenomicRanges::start(vars) %% 2 == 1
    vars[odd.pos]$baf <- 1 - vars[odd.pos]$baf
  }

  vars <- sort(vars)

  return(vars)
}


