#' plotLRR
#'
#' @description
#' Plots the raw SNP array data using karyoploteR
#'
#' @details
#' Creates a plot with the LRR values along the genome
#'
#' @usage plotLRR(karyoplot, snps, lrr.column="lrr", label="LRR", ymin=-4, ymax=2, out.of.range = "points", out.of.range.col="red", density.height=0.05, density.window=1e5,
#'                line.at.0 = TRUE, line.at.0.col="blue", r0=0, r1=1, points.cex=0.3, points.col="#333333",
#'                points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, axis.cex=1.2, ...)
#'
#' @param karyoplot the karyoplot
#' @param snps the data
#' @param lrr.column (defaults to "lrr")
#' @param label (defaults to "LRR")
#' @param ymin (defaults to -4)
#' @param ymax (defaults to 2)
#' @param out.of.range  (defaults to  "points")
#' @param out.of.range.col (defaults to "red")
#' @param density.height (defaults to 0.05)
#' @param density.window (defaults to 1e5)
#' @param line.at.0  (defaults to  TRUE)
#' @param line.at.0.col (defaults to "blue")
#' @param r0 (defaults to 0)
#' @param r1 (defaults to 1)
#' @param points.cex (defaults to 0.3)
#' @param points.col (defaults to "#333333")
#' @param points.pch (defaults to 16)
#' @param label.cex (defaults to 1.5)
#' @param label.srt (defaults to 90)
#' @param label.margin (defaults to 0.03)
#' @param axis.cex (defaults to 1.2)
#' @param ... Additional params
#'
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#'
#' pos <- floor(runif(1000, 1, 10000))
#' lrr.data <- toGRanges("chr1", pos, pos)
#' lrr.data$lrr <- rnorm(1000, mean = 0, sd = 0.3)
#'
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotLRR(kp, lrr.data)
#'
#'
#' names(mcols(lrr.data)) <- "values"
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotLRR(kp, lrr.data, lrr.column="values", r0=0, r1=0.5, points.col="red", label="Values")
#' plotLRR(kp, lrr.data, lrr.column=1, r0=0.5, r1=1, points.col="orange", label="First", label.cex=0.8, points.cex=1.4)
#'
#'
#' @export plotLRR
#'
#' @importFrom stats rnorm
#'

plotLRR <- function(karyoplot, snps, lrr.column="lrr", label="LRR", ymin=-4, ymax=2, out.of.range = "points", out.of.range.col="red", density.height=0.05, density.window=1e5,
                    line.at.0 = TRUE, line.at.0.col="blue",
                    r0=0, r1=1, points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, axis.cex=1.2, ...) {

  snps <- regioneR::toGRanges(snps)


  karyoploteR::kpAddLabels(karyoplot, r0=r0, r1=r1, labels = label, srt=label.srt, pos = 3, cex=label.cex, label.margin = label.margin, ...)
  if(ymin<0 && ymax>0) {
    tick.pos <- c(ceiling(ymin), 0, floor(ymax))
  } else {
    tick.pos <- c(ceiling(ymin), floor(ymax))
  }
  karyoploteR::kpAxis(karyoplot, r0=r0, r1=r1, ymin=ymin, ymax=ymax, tick.pos = tick.pos, cex=axis.cex, ...)


  #Plot the points
  #Identify out of range points
  below.min <- GenomicRanges::mcols(snps)[,lrr.column] < ymin
  above.max <- GenomicRanges::mcols(snps)[,lrr.column] > ymax

  #plot the "in range" points
  karyoploteR::kpPoints(karyoplot, data=snps[!(below.min | above.max)], y=GenomicRanges::mcols(snps)[,lrr.column][!(below.min | above.max)],
           ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=points.col, pch=points.pch, cex=points.cex, ...)

  #and plot the out of range points
  if(out.of.range == "points") {
    #Plot them as individual points with a bit of jitter to make them more visible
    if(any(above.max)) {
      karyoploteR::kpPoints(karyoplot, data=snps[above.max], y=ymax - stats::rnorm(length(which(above.max)), 0.02, 0.01), ymin=ymin, ymax=ymax,
               r0=r0, r1=r1, pch=points.pch, cex=points.cex, col=out.of.range.col, ...)
    }
    if(any(below.min)) {
      karyoploteR::kpPoints(karyoplot, data=snps[below.min], y=ymin - stats::rnorm(length(which(below.min)), 0.02, 0.01), ymin=ymin, ymax=ymax,
               r0=r0, r1=r1, pch=points.pch, cex=points.cex, col=out.of.range.col, ...)
    }
  } else if(out.of.range == "density") {
    #Idea: make the density plots face outward of the main plot, reducing the height of the points plot by 2*density.height
    if(any(above.max)) {
      karyoploteR::kpPlotDensity(karyoplot, data=snps[above.max], r0=r1, r1=r1-density.height*abs((r1-r0)),
                    window.size = density.window, col=out.of.range.col, border=NULL, ...)
    }
    if(any(below.min)) {
      karyoploteR::kpPlotDensity(karyoplot, data=snps[below.min], r0=r0, r1=r0+density.height*abs((r1-r0)),
                    window.size = density.window, col=out.of.range.col, border=NULL, ...)
    }
  } else {
  }

  if(line.at.0==TRUE) {
    karyoploteR::kpAbline(karyoplot, h=0, ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=line.at.0.col)
  }

  invisible(karyoplot)
}
