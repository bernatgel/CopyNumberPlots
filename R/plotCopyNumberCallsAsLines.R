#' plotCopyNumberCallsAsLines
#'
#' @description
#' Plot the segments representing the copy number calls by any algorithm
#'
#' @details
#' Plots the segments
#'
#' @usage plotCopyNumberCallsAsLines(karyoplot, cn.calls, style="line", cn.column="cn", label="", lab.cex=1, axis=TRUE, axis.cex=1, numticks=NULL, col="black", ymin=NULL, ymax=NULL, r0=0, r1=1, ...)
#'
#' @param karyoplot The karyoplot
#' @param cn.calls  The data
#' @param style (defaults to "line")
#' @param cn.column (defaults to "cn")
#' @param label (defaults to "")
#' @param lab.cex (defaults to 1)
#' @param axis (defaults to TRUE)
#' @param axis.cex (defaults to 1)
#' @param numticks (defaults to NULL)
#' @param col (defaults to "black")
#' @param ymin (defaults to NULL)
#' @param ymax (defaults to NULL)
#' @param r0 (defaults to 0)
#' @param r1 (defaults to 1)
#' @param ... additional params
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
#' plotCopyNumberCallsAsLines(kp, scnas, r0=0.2, r1=0.4, style="segments", axis=FALSE, label="CopyNumber")
#'
#' @export plotCopyNumberCallsAsLines
#'
#'



plotCopyNumberCallsAsLines <- function(karyoplot, cn.calls, style="line", cn.column="cn", label="", lab.cex=1, axis=TRUE, axis.cex=1, numticks=NULL, col="black", ymin=NULL, ymax=NULL, r0=0, r1=1, ...) {


  #If cn.calls is a list, call this same function with each single element to actually produce the plot. Use autotrack to set the appropiate r0 and r1 values.
  if(methods::is(cn.calls, "list")) {
    for(i in seq_len(length(cn.calls))) {
      if(methods::is(labels, "list")) labels <- as.character(unlist(labels))
      if(length(labels)==length(cn.calls)) { #If there are as many labels as samples, assume each label should be used for one track
        lab <- labels[i]
      } else {
        lab <- labels[1]
      }
      rr <- karyoploteR::autotrack(current.track=i, total.tracks=length(cn.calls), r0=r0, r1=r1)
      plotCopyNumberCallsAsLines(karyoplot, cn.calls[[i]], r0=rr$r0, r1=rr$r1,
                          labels = lab, lab.cex = lab.cex, style=style,
                          cn.column=cn.column, axis=axis, axis.cex=axis.cex, numticks=numticks,
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

    karyoploteR::kpLines(karyoplot, data=cn.calls, y=GenomicRanges::mcols(cn.calls)[,cn.column], ymin=ymin, ymax=ymax, col=col, ...)
  } else { #Plot as segments (without the vertical lines between them)
    karyoploteR::kpSegments(karyoplot, data=cn.calls, y0=GenomicRanges::mcols(cn.calls)[,cn.column], y1=GenomicRanges::mcols(cn.calls)[,cn.column], ymin=ymin, ymax=ymax, col=col,...)
  }

  if(axis==TRUE) {
    if(is.null(numticks)) numticks <- ymax-ymin+1
    karyoploteR::kpAxis(karyoplot, ymin=ymin, ymax=ymax, cex=axis.cex, numticks=numticks, ...)
  }
  if(!is.null(label) & nchar(label)>0) {
    karyoploteR::kpAddLabels(karyoplot, labels = label, cex=lab.cex, ...)
  }

  invisible(karyoplot)
}

