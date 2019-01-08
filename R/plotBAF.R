#' plotBAF
#'
#' @description Plot the B-allele frequency (BAF) data
#'
#' @details This function plot the B-allele frequency (BAF) values on the
#' genome. BAF values represent the frequency of one of the alleles (NOT always
#' the minor allele) in the sample and usually come from SNP-array. However,
#' it's possible to get similar values from NGS variant calling data.
#' The function plots the data points, an optional axis and an optional label identifying the data.
#' The function expects a GRanges object with the position of the SNPs with a column
#' (usually named 'baf') with the BAF values or any object valid to \code{\link[regioneR]{toGRanges}}. This object can be created with
#' \code{\link{loadSNPData}}.
#'
#' @usage plotBAF(karyoplot, snps, baf.column="baf", label="BAF", points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, add.axis=TRUE, axis.cex=1.2, r0=0, r1=1, data.panel=1, ...)
#'
#'
#' @param karyoplot  (a KaryoPlot object) The object returned by the \code{\link[karyoploteR]{plotKaryotype}} function and representing the current active plot.
#' @param snps  (a GRanges) A GRanges object with the positions of the SNPs and a column with the BAF values. Other columns are ignored.
#' @param baf.column (a character) The name of the column in \code{snps} with the BAF values. (defaults to "baf")
#' @param label (a character) The text of the label to identify the data. If NULL or NA, no label will be plotted. (defaults to "BAF")
#' @param points.cex (numeric) The size of the points. (defaults to 0.3)
#' @param points.col (a color) The color of the points (defaults to "#333333")
#' @param points.pch (numeric) The shape of the points. (defaults to 16, a filled circle)
#' @param label.cex (numeric) The size of the label (defaults to 1.5)
#' @param label.srt (numeric) The rotation of the label (defaults to 90, vertical text)
#' @param label.margin (numeric) The margin between the label and the origin of the plot in plot coordinates (the width of the plot is 1). (defaults to 0.03)
#' @param add.axis (logical) Whether to add an axis (defaults to TRUE)
#' @param axis.cex (numeric) The size of the axis labels.  (defaults to 1.2)
#' @param data.panel    (numeric) (karyoploteR parameter) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to karyoploteR functions.
#'
#'
#'
#' @return Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#'
#' pos <- sort(floor(runif(1000, 1, 10000)))
#' baf.data <- toGRanges("chr1", pos, pos)
#' baf.data$baf <- rnorm(1000, mean = 0.5, sd = 0.05)
#' baf.data$baf[1:400] <- baf.data$baf[1:400] + c(0.2, -0.2)
#'
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotBAF(kp, baf.data)
#'
#'
#' names(mcols(baf.data)) <- "values"
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotBAF(kp, baf.data, baf.column="values", r0=0, r1=0.5, points.col="red", label="Values")
#' plotBAF(kp, baf.data, baf.column=1, r0=0.5, r1=1, points.col="orange", label="First", label.cex=0.8, points.cex=1.4)
#'
#' @export plotBAF
#'
#' @importFrom GenomicRanges mcols



plotBAF <- function(karyoplot, snps, baf.column="baf", label="BAF", points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, add.axis=TRUE, axis.cex=1.2, r0=0, r1=1, data.panel=1, ...) {


  #TODO: Add  parameter checks

  snps <- regioneR::toGRanges(snps)

  if(!is.null(label) && !is.na(label) && nchar(label)>0) {
    karyoploteR::kpAddLabels(karyoplot, r0=r0, r1=r1, labels = label, srt=label.srt, pos = 3, cex=label.cex, label.margin = label.margin, data.panel=data.panel, ...)
  }

  if(add.axis==TRUE) {
    karyoploteR::kpAxis(karyoplot, r0=r0, r1=r1, ymin=0, ymax=1, cex=axis.cex, data.panel=data.panel, ...)
  }

  karyoploteR::kpPoints(karyoplot, data=snps, y=GenomicRanges::mcols(snps)[,baf.column], ymin=0, ymax=1, r0=r0, r1=r1, col=points.col, pch=points.pch, cex=points.cex, data.panel=data.panel, ...)

  invisible(karyoplot)

}
