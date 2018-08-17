#' plotSNPData
#'
#' @description
#' Plots the raw SNP array data using karyoploteR
#'
#' @details
#' Creates a plot with the LRR and BAF values along the genome
#'
#' @usage plotSNPData(snps, main="Raw Data", chromosomes="canonical", zoom=NULL, lrr.min=-4, lrr.max=2, total.height=1, bottom=0, margin=0.05, points.cex=0.3, labels.cex=1.5, main.cex=2, axis.cex=1.2, chr.cex=1.5)
#'
#' @param snps The SNP array data
#' @param main (defaults to "Raw Data")
#' @param chromosomes (defaults to "canonical")
#' @param zoom (defaults to NULL)
#' @param lrr.min (defaults to -4)
#' @param lrr.max (defaults to 2)
#' @param total.height (defaults to 1)
#' @param bottom (defaults to 0)
#' @param margin (defaults to 0.05)
#' @param points.cex (defaults to 0.3)
#' @param labels.cex (defaults to 1.5)
#' @param main.cex (defaults to 2)
#' @param axis.cex (defaults to 1.2)
#' @param chr.cex (defaults to 1.5)
#'
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#'
#'
#'
#' @export plotBAF
#'
#'


plotBAF <- function(karyoplot, snps, label="BAF", r0=0, r1=1, points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, axis.cex=1.2) {

  snps <- regioneR::toGRanges(snps)
  if(!("baf" %in% names(mcols(snps)))) stop("The snps object need a column named 'baf'.")

  kpAddLabels(karyoplot, r0=r0, r1=r1, labels = label, srt=label.srt, pos = 3, cex=label.cex, label.margin = label.margin)
  kpAxis(karyoplot, r0=r0, r1=r1, ymin=0, ymax=1, cex=axis.cex)
  kpPoints(karyoplot, data=snps, y=snps$baf, ymin=0, ymax=1, r0=r0, r1=r1, col=points.col, pch=points.pch, cex=points.cex)

  invisible(karyoplot)

}
