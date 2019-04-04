#' plotCopyNumberCallsAsLines
#'
#' @description
#' Plot the segments representing the copy number calls by any algorithm
#'
#' @details
#' Plots the segments
#'
#' @usage plotCopyNumberCallsAsLines(karyoplot, cn.calls, style="line", cn.column="cn", labels=NULL, label.cex=1, add.axis=TRUE, axis.cex=1, numticks=NULL, col="black", ymin=NULL, ymax=NULL, r0=0, r1=1, ...)
#'
#' @param karyoplot (a KaryoPlot object) The object returned by the \code{\link[karyoploteR]{plotKaryotype}} function and representing the current active plot.
#' @param cn.calls  (a GRanges, a list of GRanges or a GRangesList) An object with the positions of the CN calls and a column with the CN values. Other columns are ignored. If it's a list of GRanges with different samples, all samples will be plotted, splitting the total plot space between them.
#' @param style (character) The style in which the lines can be plot. It colud be as lines or segments using the kplines or kpsegments. (defaults to "line")
#' @param cn.column (integer or character vector) The name or number of the column with CN information.(defaults to "cn")
#' @param labels (character vector or list) labels (a character) The text of the label to identify the data. If NA, no label will be plotted. If NULL, if snps is a single sample GRanges it will default to "BAF", if it's a list of samples it will default to the names in the list or consecutive numbers if names(snps) is NULL. (defaults to NULL)
#' @param label.cex (numeric) The size of the label (defaults to 1)
#' @param axis (logical) Whether to plot an axis (defaults to TRUE)
#' @param axis.cex (numeric) The size of the axis labels.(defaults to 1)
#' @param numticks (defaults to NULL)
#' @param col (color) The color of the lines (defaults to "black")
#' @param ymin (numeric) (karyoploteR parameter) The minimum value of y to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to -4)
#' @param ymax (numeric) (karyoploteR parameter) (numeric) The maximum value of y to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to 2)
#' @param r0  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to 0)
#' @param r1  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to 1)
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
#' set.seed(1234)
#' gg <- filterChromosomes(getGenome("hg19"))
#'
#' scnas <- createRandomRegions(40, 10e6, 10e6, non.overlapping=TRUE)
#' scnas$cn <- floor(runif(40, 0, 4))
#' scnas$loh <- ifelse(scnas$cn<2, 1, 0)
#' normal.regs <- subtractRegions(gg, scnas)
#' normal.regs$cn <- 2
#' normal.regs$loh <- 0
#' scnas <- sort(c(scnas, normal.regs))
#'
#' kp <- plotKaryotype("hg19")
#' plotCopyNumberCallsAsLines(kp, scnas)
#'
#' kp <- plotKaryotype("hg19", chromosomes="chr3")
#' plotCopyNumberCalls(kp, scnas, r0=0, r1=0.2)
#' plotCopyNumberCallsAsLines(kp, scnas, r0=0.25, r1=0.55)
#' plotCopyNumberCallsAsLines(kp, scnas, r0=0.65, r1=0.95, style="segments", col="red", ymin=0, ymax=4, numticks=3)
#'
#' kp <- plotKaryotype("hg19")
#' plotCopyNumberCallsAsLines(kp, scnas, r0=0.2, r1=0.4, style="segments", add.axis=FALSE, labels="CopyNumber")
#'
#' @export plotCopyNumberCallsAsLines
#'
#'



plotCopyNumberCallsAsLines <- function(karyoplot, cn.calls, style="line", cn.column="cn", labels=NULL, label.cex=1, add.axis=TRUE, axis.cex=1, numticks=NULL, col="black", ymin=NULL, ymax=NULL, r0=0, r1=1, ...) {

  #If cn.calls is a list, call this same function with each single element to actually produce the plot. Use autotrack to set the appropiate r0 and r1 values.
  if(methods::is(cn.calls, "list") || methods::is(cn.calls, "CompressedGRangesList")) {
    if(is.null(labels)) labels <- ifelse(is.null(names(cn.calls)), seq_len(length(cn.calls)), names(cn.calls))
    if(methods::is(labels, "list")) labels <- as.character(unlist(labels)) #If labels are a list, make them a vector
    for(i in seq_len(length(cn.calls))) {
      #If there are as many labels as samples, assume each label should be used for one track, else, use the first one
      lab <- ifelse(length(labels)==length(cn.calls), labels[i], labels[1])
      rr <- karyoploteR::autotrack(current.track=i, total.tracks=length(cn.calls), margin = track.margin, r0=r0, r1=r1)
      plotCopyNumberCallsAsLines(karyoplot, cn.calls[[i]], r0=rr$r0, r1=rr$r1,
                                 labels = lab, label.cex = label.cex, style=style,
                                 cn.column=cn.column, add.axis=add.axis, axis.cex=axis.cex, numticks=numticks,
                                 col=col, ymin=ymin, ymax=ymax, ...)
    }
    return(invisible(karyoplot))
  }
  

  style <- match.arg(style, choices = c("line", "segments"))

  if(is.null(ymin)) ymin <- min(GenomicRanges::mcols(cn.calls)[,cn.column])
  if(is.null(ymax)) ymax <- max(GenomicRanges::mcols(cn.calls)[,cn.column])

  if(style=="line") {
    starts <- cn.calls
    GenomicRanges::end(starts) <- GenomicRanges::start(starts)
    ends <- cn.calls
    GenomicRanges::start(ends) <- GenomicRanges::end(ends)
    cn.calls <- sort(c(starts, ends))

    cn.calls <- sort(cn.calls)

    karyoploteR::kpLines(karyoplot, data=cn.calls, y=GenomicRanges::mcols(cn.calls)[,cn.column], ymin=ymin, ymax=ymax, col=col, r0=r0, r1=r1, ...)
  } else { #Plot as segments (without the vertical lines between them)
    karyoploteR::kpSegments(karyoplot, data=cn.calls, y0=GenomicRanges::mcols(cn.calls)[,cn.column], y1=GenomicRanges::mcols(cn.calls)[,cn.column], ymin=ymin, ymax=ymax, col=col, r0=r0, r1=r1, ...)
  }

  if(add.axis==TRUE) {
    if(is.null(numticks)) numticks <- ymax-ymin+1
    karyoploteR::kpAxis(karyoplot, ymin=ymin, ymax=ymax, cex=axis.cex, numticks=numticks, r0=r0, r1=r1, ...)
  }
  
  if(methods::is(labels, "list")) labels <- as.character(unlist(labels)) #If labels are a list, make them a vector
  
  if(is.null(labels)) labels <- "Segments"
  
  if(!is.na(labels) && nchar(labels)>0) {
    karyoploteR::kpAddLabels(karyoplot, labels = labels, cex=label.cex, r0=r0, r1=r1, ...)
  }
  

  invisible(karyoplot)
}

