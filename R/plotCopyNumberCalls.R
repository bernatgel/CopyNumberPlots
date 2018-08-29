#' plotCopyNumberCalls
#'
#' @description
#' Plot the segments representing the copy number calls by any algorithm
#'
#' @details
#' Plots the segments
#'
#' @usage plotCopyNumberCalls()
#'
#' @param karyoplot
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
#' @export plotCopyNumberCalls
#'
#' @import karyoploteR
#' @importFrom IRanges subsetByOverlaps
#'





plotCopyNumberCalls <- function(karyoplot, cn.calls, cn.colors=NULL, loh.color="#1E90FF", labels="", lab.cex=1, lab2.cex=NULL, loh.height=0.3, r0=0, r1=1) {

  if(is.null(lab2.cex)) lab2.cex <- lab.cex

  segment.colors <- getCopyNumberColors(cn.colors)

  #Explicitly filter the segments, since it will use the wrong colors otherwise
  #Check: is it still neede? The bug should be fixed
  plot.region <- karyoplot$plot.region
  segments <- IRanges::subsetByOverlaps(cn.calls, plot.region)

  seg.cols <- segment.colors[as.character(cn.calls$cn)]


  kpRect(kp, data=cn.calls, y0=loh.height, y1=1, col=seg.cols, r0=r0, r1=r1, border=NA)
  kpRect(kp, data=cn.calls[cn.calls$loh==1], y0=0, y1=loh.height, r0=r0, r1=r1, col=loh.color, border=NA)
  if(!is.null(labels) && !is.na(labels) && all(is.character(labels) && length(labels)>0)) {
    if(length(labels)==1) {
      kpAddLabels(kp, labels = labels[1], r0=r0, r1=r1, cex=lab.cex)
    } else { #If there are two labels, use one for the CN track an another for the LOH track
      kpAddLabels(kp, labels = labels[1], r0=r0+(r1-r0)*loh.height, r1=r1, cex=lab.cex)
      kpAddLabels(kp, labels = labels[2], r0=r0, r1=r0+(r1-r0)*loh.height, cex=lab2.cex)
    }
  }
  invisible(karyoplot)
}

