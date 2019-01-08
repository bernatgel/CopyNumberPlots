#' plotCopyNumberSummary
#'
#' @description
#' Plot a summary of copy number status over a number of samples using a
#' histogram-like representation.
#'
#' @details
#' 2 types of plots. Convergent or divergent barplots
#'
#' @usage plotCopyNumberSummary(karyoplot, cn.calls, direction="in",  gain.color=NULL, normal.color=NULL, loss.color=NULL, add.grid=FALSE, grid.color="white", label="", lab.cex=1, lab.srt=0, pos=2, r0=0, r1=1, ...)
#'
#' @param karyoplot the karyoplot object
#' @param cn.calls (list) a LIST of cn.calls objects!
#' @param direction (defaults to "in")
#' @param gain.color (defaults to NULL)
#' @param normal.color (defaults to NULL)
#' @param loss.color (defaults to NULL)
#' @param add.grid (defaults to FALSE)
#' @param grid.color (defaults to "white")
#' @param label (defaults to "")
#' @param lab.cex (defaults to 1)
#' @param lab.srt (defaults to 0)
#' @param pos (defaults to 2)
#' @param r0 (defaults to 0)
#' @param r1 (defaults to 1)
#' @param ... additional params
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#'
#'
#' gg <- filterChromosomes(getGenome("hg19"))
#'
#' all.scnas <- list()
#' for(i in seq_len(10)) {
#'   scnas <- createRandomRegions(40, 10e6, 10e6)
#'   scnas$cn <- floor(runif(40, 0, 4))
#'   scnas$loh <- ifelse(scnas$cn<2, 1, 0)
#'   normal.regs <- subtractRegions(gg, scnas)
#'   normal.regs$cn <- 2
#'   scnas <- sort(c(scnas, normal.regs))
#'   all.scnas[[paste0("samp", i)]] <- scnas
#' }
#'
#' kp <- plotKaryotype("hg19", plot.type=4)
#' plotCopyNumberSummary(kp, all.scnas)
#'
#' kp <- plotKaryotype("hg19", plot.type=4)
#' plotCopyNumberSummary(kp, all.scnas, gain.col="blue", loss.col="red", label="Copy Number", r0=0, r1=0.3, lab.srt=90, pos=3)
#' plotCopyNumberSummary(kp, all.scnas, direction="out", add.grid=TRUE, r0=0.35, r1=0.65)
#' plotCopyNumberSummary(kp, all.scnas, direction="out", r0=0.7, r1=1)
#' kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=1, tick.pos = c(0,5,10))
#' kpAxis(kp, ymin=0, ymax=10, r0=0.85, r1=0.7, tick.pos = c(5,10))
#'
#'
#' @export plotCopyNumberSummary
#'
#' @importFrom GenomicRanges GRangesList


plotCopyNumberSummary <- function(karyoplot, cn.calls, direction="in",  gain.color=NULL, normal.color=NULL, loss.color=NULL, add.grid=FALSE, grid.color="white", label="", lab.cex=1, lab.srt=0, pos=2, r0=0, r1=1, ...) {

  if(!is.list(cn.calls)) stop("cn.calls must be a list of GRanges with copy number information")

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

  if(nchar(label)>0) {
    karyoploteR::kpAddLabels(karyoplot, labels = label, r0=r0, r1=r1, cex=lab.cex, srt=lab.srt, pos=pos, ...)
  }

  invisible(karyoplot)
}

