#' plotLRR
#'
#' @description
#' Plots the raw SNP array data using karyoploteR
#'
#' @details
#' Creates a plot with the LRR values along the genome
#'
#' @usage plotLRR(karyoplot, snps, lrr.column="lrr", labels="LRR", ymin=-4, ymax=2, out.of.range = "points", out.of.range.col="red", density.height=0.05, density.window=1e5, line.at.0 = TRUE, line.at.0.col="blue", r0=0, r1=1, points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, add.axis=TRUE, axis.cex=1.2, track.margin=0.1, data.panel=1, ...)
#'
#' @param karyoplot (a KaryoPlot object) The object returned by the \code{\link[karyoploteR]{plotKaryotype}} function and representing the current active plot.
#' @param snps (a GRanges, a list of GRanges or a GRangesList) An object with the positions of the SNPs and a column with the BAF values. Other columns are ignored. If it's a list of GRanges with different samples, all samples will be plotted, splitting the total plot space between them.
#' @param lrr.column (number or character) The name or number of the column with LRR information. (defaults to "lrr")
#' @param labels (character) The text of the label to identify the data. If NA, no label will be plotted. If NULL, if snps is a single sample GRanges it will default to "BAF", if it's a list of samples it will default to the names in the list or consecutive numbers if names(snps) is NULL.(defaults to "LRR")
#' @param ymin (numeric) (karyoploteR parameter) The minimum value of y to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to -4)
#' @param ymax (numeric) (karyoploteR parameter) (numeric) The maximum value of y to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to 2)
#' @param out.of.range (a character) Either to plot "points" or "density" (defaults to  "points")
#' @param out.of.range.col (defaults to "red")
#' @param density.height (defaults to 0.05)
#' @param density.window (defaults to 1e5)
#' @param line.at.0  (logical) Whether to plot an horizontal line at 0. (defaults to  TRUE)
#' @param line.at.0.col (color) The color of the horizontal line plotted at 0. (defaults to "blue")
#' @param r0  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to 0)
#' @param r1  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to 1)
#' @param points.cex (numeric) The size of the points. (defaults to 0.3)
#' @param points.col (a color) The color of the points(defaults to "#333333")
#' @param points.pch (numeric) The shape of the points. (defaults to 16, a filled circle)
#' @param label.cex (numeric) The size of the label. (defaults to 1.5)
#' @param label.srt (numeric) The rotation of the label. (defaults to 90)
#' @param label.margin (numeric) The margin between the label and the origin of the plot in plot coordinates (the width of the plot is 1). (defaults to 0.03)
#' @param add.axis (logical) Whether to add an axis (defaults to TRUE)
#' @param axis.cex (numeric) The size of the axis labels.  (defaults to 1.2)
#' @param track.margin (numeric) (defaults to 0.1)
#' @param data.panel (numeric) (karyoploteR parameter) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to karyoploteR functions.
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
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotLRR(kp, lrr.data, lrr.column=1, points.col="orange", labels="First", label.cex=0.8, points.cex=1.4)
#'
#'
#' @export plotLRR
#'
#' @importFrom stats rnorm
#'

plotLRR <- function(karyoplot, snps, lrr.column="lrr", labels="LRR", ymin=-4, ymax=2, out.of.range = "points", out.of.range.col="red", density.height=0.05, density.window=1e5,
                    line.at.0 = TRUE, line.at.0.col="blue",
                    r0=0, r1=1, points.cex=0.3, points.col="#333333", points.pch=16, label.cex=1.5, label.srt=90, label.margin=0.03, add.axis=TRUE, axis.cex=1.2, track.margin=0.1, data.panel=1, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop("karyoplot must be a KaryoPlot object")
  out.of.range <- match.arg(out.of.range, c("points", "density"))

  #If cn.calls is a list, call this same function with each single element to actually produce the plot. Use autotrack to set the appropiate r0 and r1 values.
  if(methods::is(snps, "list") || methods::is(snps, "CompressedGRangesList")) {
    if(is.null(labels)) labels <- ifelse(is.null(names(snps)), seq_len(length(snps)), names(snps))
    for(i in seq_len(length(snps))) {
      #If there are as many labels as samples, assume each label should be used for one track, else, use the first one
      lab <- ifelse(length(labels)==length(snps), labels[i], labels[1])

      rr <- karyoploteR::autotrack(current.track=i, total.tracks=length(snps), margin = track.margin, r0=r0, r1=r1)
      plotLRR(karyoplot, snps[[i]], r0=rr$r0, r1=rr$r1, labels=lab,
              lrr.column=lrr.column, ymin=ymin, ymax=ymax,
              out.of.range=out.of.range, out.of.range.col=out.of.range.col,
              density.height=density.height, density.window=density.window,
              line.at.0 = line.at.0, line.at.0.col=line.at.0.col,
              points.cex=points.cex, points.col=points.col, points.pch=points.pch,
              label.cex=label.cex, label.srt=label.srt, label.margin=label.margin,
              add.axis=add.axis, axis.cex=axis.cex, track.margin = track.margin,
              data.panel=data.panel, ...)
    }
    return(invisible(karyoplot))
  }





  snps <- regioneR::toGRanges(snps)

  if(lrr.column!="lrr") names(GenomicRanges::mcols(snps))[which(names(GenomicRanges::mcols(snps))==lrr.column)] <- "lrr"

  snps <- removeNAs(snps, lrr.na = TRUE, baf.na = FALSE, id.na = FALSE, verbose = FALSE)


  if(is.null(labels)) labels <- "LRR"

  if(!is.na(labels) && nchar(labels)>0) {
    karyoploteR::kpAddLabels(karyoplot, r0=r0, r1=r1, labels = labels, srt=label.srt, pos = 3, cex=label.cex, label.margin = label.margin, data.panel=data.panel, ...)
  }


  if(add.axis==TRUE) {
    if(ymin<0 && ymax>0) {
      tick.pos <- c(ceiling(ymin), 0, floor(ymax))
    } else {
      tick.pos <- c(ceiling(ymin), floor(ymax))
    }
    karyoploteR::kpAxis(karyoplot, r0=r0, r1=r1, ymin=ymin, ymax=ymax, tick.pos = tick.pos, cex=axis.cex, ...)
  }

  #Plot the points
  #Identify out of range points
  below.min <- GenomicRanges::mcols(snps)[,lrr.column] < ymin
  above.max <- GenomicRanges::mcols(snps)[,lrr.column] > ymax

  #plot the "in range" points
  if(any(!(below.min | above.max))) {
    karyoploteR::kpPoints(karyoplot, data=snps[!(below.min | above.max)], y=GenomicRanges::mcols(snps)[,lrr.column][!(below.min | above.max)],
           ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=points.col, pch=points.pch, cex=points.cex, ...)
  }

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
  } 

  if(line.at.0==TRUE) {
    karyoploteR::kpAbline(karyoplot, h=0, ymin=ymin, ymax=ymax, r0=r0, r1=r1, col=line.at.0.col)
  }

  invisible(karyoplot)
}
