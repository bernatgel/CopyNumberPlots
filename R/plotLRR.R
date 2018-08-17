#' plotLRR
#'
#' @description
#' Plots the raw SNP array data using karyoploteR
#'
#' @details
#' Creates a plot with the LRR values along the genome
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
#' @export plotLRR
#'
#'


plotLRR <- function(karyoplot, snps, label="LRR", ymin=-4, ymax=2, out.of.range = "points", out.of.range.col="red", density.height=0.05, density.window=1e5,
                    line.at.0 = TRUE, line.at.0.col="blue",
                    r0=0, r1=1, points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, axis.cex=1.2, ...) {

  snps <- regioneR::toGRanges(snps)
  if(!("lrr" %in% names(mcols(snps)))) stop("The snps object need a column named 'lrr'.")

  kpAddLabels(karyoplot, r0=r0, r1=r1, labels = label, srt=label.srt, pos = 3, cex=label.cex, label.margin = label.margin, ...)
  if(ymin<0 && ymax>0) {
    tick.pos <- c(ceiling(ymin), 0, floor(ymax))
  } else {
    tick.pos <- c(ceiling(ymin), floor(ymax))
  }
  kpAxis(karyoplot, r0=r0, r1=r1, ymin=ymin, ymax=ymax, tick.pos = tick.pos, cex=axis.cex, ...)


  #Plot the points
  #Identify out of range points
  below.min <- snps$lrr < ymin
  above.max <- snps$lrr > ymax

  #plot the "in range" points
  kpPoints(karyoplot, data=snps[!(below.min | above.max)], y=snps$lrr[!(below.min | above.max)],
           ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=points.col, pch=points.pch, cex=points.cex, ...)

  #and plot the out of range points
  if(out.of.range == "points") {
    #Plot them as individual points with a bit of jitter to make them more visible
    if(any(above.max)) {
      kpPoints(karyoplot, data=snps[above.max], y=ymax - rnorm(length(which(above.max)), 0.02, 0.01), ymin=ymin, ymax=ymax,
               r0=r0, r1=r1, pch=points.pch, cex=points.cex, col=out.of.range.col, ...)
    }
    if(any(below.min)) {
      kpPoints(karyoplot, data=snps[below.min], y=ymin - rnorm(length(which(below.min)), 0.02, 0.01), ymin=ymin, ymax=ymax,
               r0=r0, r1=r1, pch=points.pch, cex=points.cex, col=out.of.range.col, ...)
    }
  } else if(out.of.range == "density") {
    #Idea: make the density plots face outward of the main plot, reducing the height of the points plot by 2*density.height
    if(any(above.max)) {
      kpPlotDensity(karyoplot, data=snps[above.max], r0=r1, r1=r1-density.height*abs((r1-r0)),
                    window.size = density.window, col=out.of.range.col, border=NULL, ...)
    }
    if(any(below.min)) {
      kpPlotDensity(karyoplot, data=snps[below.min], r0=r0, r1=r0+density.height*abs((r1-r0)),
                    window.size = density.window, col=out.of.range.col, border=NULL, ...)
    }
  } else {
    #If not points or density, assume "hide" and don't plot them
  }

  if(line.at.0==TRUE) {
    kpAbline(karyoplot, h=0, ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=line.at.0.col)
  }

  invisible(karyoplot)
}
