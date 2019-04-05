#' plotCopyNumberSummary
#'
#' @description
#' Plot a summary of copy number status over a number of samples using a
#' histogram-like representation.
#'
#' @details
#' 2 types of plots. Convergent or divergent barplots
#'
#' @usage plotCopyNumberSummary(karyoplot, cn.calls, direction="in",  gain.color=NULL, normal.color=NULL, loss.color=NULL, add.grid=FALSE, grid.color="white", labels=NULL, label.cex=1, label.srt=0, pos=2, r0=0, r1=1, ...)
#'
#' @param karyoplot (a KaryoPlot object) The object returned by the \code{\link[karyoploteR]{plotKaryotype}} function and representing the current active plot.
#' @param cn.calls (a list of GRanges or a GRangesList) A list of GRanges or a GRangesList containing the GRanges objects with cn.column value
#' @param direction The direction to which the coverage plot point, either "in" for inward or "out" for outward. (defaults to "in")
#' @param gain.color (color) The color assigned to gains (defaults to NULL)
#' @param normal.color (color) The color assigned to normal ploidy(defaults to NULL)
#' @param loss.color (colors) The color assigned to losses(defaults to NULL)
#' @param add.grid (logical) Whether to add lines as a grid in the plot (defaults to FALSE)
#' @param grid.color (color) The color of the grid (defaults to "white")
#' @param labels (a character) The text of the label to identify the data. If NA, no label will be plotted. If NULL, if snps is a single sample GRanges it will default to "CN", if it's a list of samples it will default to the names in the list or consecutive numbers if names(snps) is NULL. (defaults to NULL)
#' @param label.cex (numeric) The size of the label (defaults to 1.5)
#' @param label.srt (numeric) The rotation of the label (defaults to 90, vertical text)
#' @param pos The position of the label (defaults to 2)
#' @param r0  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1  (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to karyoploteR functions.
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#' 
#' all.scnas <- list(
#'              loadCopyNumberCalls(system.file("extdata", "S0001.ASCAT.segments.txt", package = "CopyNumberPlots", mustWork = TRUE)),
#'              loadCopyNumberCalls(system.file("extdata", "S0002.ASCAT.segments.txt", package = "CopyNumberPlots", mustWork = TRUE))
#'              )
#'              
#' kp <- plotKaryotype("hg19", plot.type=4)
#' plotCopyNumberSummary(kp, all.scnas)
#' 
#' kp <- plotKaryotype("hg19", plot.type=4)
#' plotCopyNumberSummary(kp, all.scnas, gain.col="blue", loss.col="red", labels="Copy Number", r0=0, r1=0.3, label.srt=90, pos=3)
#' plotCopyNumberSummary(kp, all.scnas, direction="out", add.grid=TRUE, r0=0.35, r1=0.65)
#' plotCopyNumberSummary(kp, all.scnas, direction="out", r0=0.7, r1=1)
#' kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=1, tick.pos = c(0,5,10))
#' kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=0.7, tick.pos = c(5,10))
#' 
#' 
#' #NOT RUN (time constraints in Bioconductor checks): 
#' #  Use a random example with more realistic results
#' # gg <- filterChromosomes(getGenome("hg19"))
#'
#' # all.scnas <- list()
#' # for(i in seq_len(10)) {
#' #   scnas <- createRandomRegions(40, 10e6, 10e6)
#' #   scnas$cn <- floor(runif(40, 0, 4))
#' #   scnas$loh <- ifelse(scnas$cn<2, 1, 0)
#' #   normal.regs <- subtractRegions(gg, scnas)
#' #   normal.regs$cn <- 2
#' #   scnas <- sort(c(scnas, normal.regs))
#' #   all.scnas[[paste0("samp", i)]] <- scnas
#' # }
#'
#' # kp <- plotKaryotype("hg19", plot.type=4)
#' # plotCopyNumberSummary(kp, all.scnas)
#'
#' # kp <- plotKaryotype("hg19", plot.type=4)
#' # plotCopyNumberSummary(kp, all.scnas, gain.col="blue", loss.col="red", labels="Copy Number", r0=0, r1=0.3, label.srt=90, pos=3)
#' # plotCopyNumberSummary(kp, all.scnas, direction="out", add.grid=TRUE, r0=0.35, r1=0.65)
#' # plotCopyNumberSummary(kp, all.scnas, direction="out", r0=0.7, r1=1)
#' # kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=1, tick.pos = c(0,5,10))
#' # kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=0.7, tick.pos = c(5,10))
#'
#'
#' @export plotCopyNumberSummary
#'
#' @importFrom GenomicRanges GRangesList


plotCopyNumberSummary <- function(karyoplot, cn.calls, direction="in",  gain.color=NULL, normal.color=NULL, loss.color=NULL, add.grid=FALSE, grid.color="white", labels=NULL, label.cex=1, label.srt=0, pos=2, r0=0, r1=1, ...) {

  if(is.null(labels)) labels <- "CN"

  if(!is.na(labels) && nchar(labels)>0) {
    karyoploteR::kpAddLabels(karyoplot, labels = labels, r0=r0, r1=r1, cex=label.cex, srt=label.srt, pos=pos, ...)
  }

  direction <- match.arg(direction, c("in", "out"))

  if(any(is.null(gain.color), is.null(loss.color), is.null(normal.color))) {
    segment.colors <- getCopyNumberColors()
    if(is.null(gain.color)) gain.color <- segment.colors["3"]
    if(is.null(normal.color)) normal.color <- segment.colors["2"]
    if(is.null(loss.color)) loss.color <- segment.colors["1"]
  }


  all.cn <-  unlist(GenomicRanges::GRangesList(cn.calls))
  if(direction=="in") {
    karyoploteR::kpDataBackground(karyoplot, r0=r0, r1=r1, col=normal.color)
    karyoploteR::kpPlotCoverage(karyoplot, all.cn[all.cn$cn>2], col=gain.color,  r0=r1, r1=r0, ymax=length(cn.calls), ...)
    karyoploteR::kpPlotCoverage(karyoplot, all.cn[all.cn$cn<2], col=loss.color, r0=r0, r1=r1, ymax=length(cn.calls), ...)
    if(add.grid==TRUE) {
      karyoploteR::kpAbline(karyoplot, h=0:length(cn.calls), col=grid.color, r0=r0, r1=r1, ymax=length(cn.calls), ...)
    }
  } else {
    mid <- r0+(r1-r0)/2
    karyoploteR::kpPlotCoverage(karyoplot, all.cn[all.cn$cn>2], col=gain.color,  r0=mid, r1=r1, ymax=length(cn.calls), ...)
    karyoploteR::kpPlotCoverage(karyoplot, all.cn[all.cn$cn<2], col=loss.color, r0=mid, r1=r0, ymax=length(cn.calls), ...)
    if(add.grid==TRUE) {
      karyoploteR::kpAbline(karyoplot, h=0:length(cn.calls), col=grid.color, r0=mid, r1=r1, ymax=length(cn.calls), ...)
      karyoploteR::kpAbline(karyoplot, h=0:length(cn.calls), col=grid.color, r0=mid, r1=r0, ymax=length(cn.calls), ...)
    }
  }


  invisible(karyoplot)
}

