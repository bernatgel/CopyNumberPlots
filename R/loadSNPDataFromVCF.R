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
#' @usage loadSNPDataFromVCF(vcf.file, regions=NULL, genome="hg19", mirror.baf=TRUE, verbose=TRUE)
#'
#' @param vcf.file The name of the file with the data
#' @param regions (defaults to NULL)
#' @param genome (a character)The name of the genome (defaults to "hg19")
#' @param mirror.baf Flip the baf of about half the snps (the ones in odd positions in the genome) to achieve a mirror-like effect as in SNP arrays (defaults to TRUE)
#' @param verbose Wether information messages should be generated. (defaults to TRUE)
#'
#' @return
#' A GRanges object with a range per SNP
#'
#' @examples
#'
#' vcf.file <- system.file("extdata", "example.vcf.gz", package = "CopyNumberPlots", mustWork = TRUE)
#' snps <- loadSNPDataFromVCF(vcf.file)
#'
#' kp <- plotKaryotype(plot.type = 4)
#' plotBAF(kp, snps = snps, labels = names(snps))
#'
#' @export loadSNPDataFromVCF
#' @importFrom VariantAnnotation readVcf ScanVcfParam scanVcfHeader geno
#' @importFrom Rsamtools TabixFile
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges start

#Read a VCF file and extract coverage (LRR) and frequency (BAF) information from it
loadSNPDataFromVCF <- function(vcf.file, regions=NULL, genome="hg19", mirror.baf=TRUE, verbose=TRUE) {
  if(verbose==TRUE) message("Scanning file ", vcf.file, "...")

  if(!is.null(regions)) regions <- tryCatch(regioneR::toGRanges(regions), error = function(e){stop("regions must be in a valid format accepted by toGRanges.\n ", e)})

  vcf.header <- VariantAnnotation::scanVcfHeader(vcf.file)

  mode <- NULL
  if("AD" %in% row.names(VariantAnnotation::geno(vcf.header))) mode <- "AD"
  if(is.null(mode)) stop("The VCF file does not have the AD field in genotype. BAF/LRR computation from FREQ and DP still not implemented.")

  #TODO: Should we accept a GRanges to scan only specific regions of the file? (for zooming, etc...)
  if(is.null(regions)) {
    vars <- VariantAnnotation::readVcf(file=Rsamtools::TabixFile(vcf.file), genome = "hg19", param = VariantAnnotation::ScanVcfParam(info=NA, geno = "AD"))
  } else {
            #TODO: Add a which=regions
  }

  #TODO: Remove the indels? Or at least make sure the frequency is correct.
  #Yes, filter out INDELS



  samples <- colnames(vars)
  res <- list()
  for(s in samples) {
    v <- vars[,s]
    if(mode == "AD") {
      #Compute the freq
      ad <- VariantAnnotation::geno(v)$AD
      ad.ref <- unlist(lapply(ad, "[", 1))
      ad.alt <- unlist(lapply(ad, "[", 2))

      #BAF
      baf <- ad.alt/(ad.alt+ad.ref)
      baf[is.nan(baf)] <- NA

      #LRR
      lrr <- ad.alt+ad.ref+1
      lrr <- log(lrr/mean(lrr))
    }

    #flip the frequency of every other SNP to achieve the same mirroring effect
    #we see in SNP-arrays
    if(mirror.baf==TRUE) baf[seq_len(length(baf))%% 2 == 1] <- 1-baf[seq_len(length(baf))%% 2 == 1]

    #Build the GRanges
    v <- SummarizedExperiment::rowRanges(v)
    GenomicRanges::mcols(v) <- data.frame(baf=baf, lrr=lrr)
    res[[s]] <- v
  }



  return(res)
}


