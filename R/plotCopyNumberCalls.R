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
#' gg <- filterChromosomes(getGenome("hg19"))
#'
#' scnas <- createRandomRegions(40, 10e6, 10e6)
#' scnas$cn <- floor(runif(40, 0, 4))
#' normal.regs <- subtractRegions(gg, scnas)
#' normal.regs$cn <- 2
#' scnas <- sort(c(scnas, normal.regs))
#' scnas$loh <- ifelse(scnas$cn<2, TRUE, FALSE)
#'
#' kp <- plotKaryotype("hg19")
#' plotCopyNumberCalls(kp, scnas)
#'
#' kp <- plotKaryotype("hg19")
#' plotCopyNumberCalls(kp, scnas, cn.colors="red_blue")
#'
#'
#' @export plotCopyNumberCalls
#'
#' @import karyoploteR
#' @importFrom IRanges subsetByOverlaps
#'





plotCopyNumberCalls <- function(karyoplot, cn.calls, cn.colors=NULL, loh.color="#1E90FF", labels="", lab.cex=1, lab2.cex=NULL, loh.height=0.3, r0=0, r1=1, ...) {


  if(methods::is(cn.calls, "list")) {
    for(i in seq_len(length(cn.calls))) {
      if(is.list(labels) && length(labels)==length(cn.calls)) { #If labels are a list, assume each label should be used for one track
        lab <- labels[[i]]
      } else {
        lab <- labels
      }
      rr <- karyoploteR::processAutotrack(r0=r0, r1=r1, autotrack = list(i, length(cn.calls)))
      plotCopyNumberCalls(karyoplot, cn.calls[[i]], r0=rr["r0"], r1=rr["r1"],
                          cn.colors = cn.colors, loh.color = loh.color,
                          labels = lab, lab.cex = lab.cex, lab2.cex = lab2.cex,
                          loh.height = loh.height, ...)
    }
    return(invisible(karyoplot))
  }

  if(is.null(lab2.cex)) lab2.cex <- lab.cex

  segment.colors <- getCopyNumberColors(cn.colors)

  #Explicitly filter the segments, since it will use the wrong colors otherwise
  #Check: is it still needed? The bug should be fixed
  plot.region <- karyoplot$plot.region
  segments <- IRanges::subsetByOverlaps(cn.calls, plot.region)

  seg.cols <- segment.colors[as.character(cn.calls$cn)]


  kpRect(karyoplot, data=cn.calls, y0=loh.height, y1=1, col=seg.cols, r0=r0, r1=r1, border=NA, ...)
  if("loh" %in% names(mcols(cn.calls))) {
    #convert loh=NA to no LOH
    cn.calls$loh[is.na(cn.calls$loh)] <- FALSE
    kpRect(karyoplot, data=cn.calls[cn.calls$loh==TRUE], y0=0, y1=loh.height, r0=r0, r1=r1, col=loh.color, border=NA, ...)
  }
  if(!is.null(labels) && !is.na(labels) && all(is.character(labels) && length(labels)>0)) {
    if(length(labels)==1) {
      kpAddLabels(karyoplot, labels = labels[1], r0=r0, r1=r1, cex=lab.cex, ...)
    } else { #If there are two labels, use the first one for the CN track an the second one for the LOH track
      kpAddLabels(karyoplot, labels = labels[1], r0=r0+(r1-r0)*loh.height, r1=r1, cex=lab.cex, ...)
      kpAddLabels(karyoplot, labels = labels[2], r0=r0, r1=r0+(r1-r0)*loh.height, cex=lab2.cex, ...)
    }
  }
  invisible(karyoplot)
}

