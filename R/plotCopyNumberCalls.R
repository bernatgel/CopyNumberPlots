#' plotCopyNumberCalls
#'
#' @description
#' Plot the segments representing the copy number calls by any algorithm
#'
#' @details
#' Plots the copy number calls as colored rectangles in the karyoplot. Each
#' copy number status has its own color (including 2n, normal genome) and
#' these colors may be specified by cn.colors. If LOH data is available
#' it will be plotted below the copy number data, with a colored rectangle
#' for every LOH region. The input is simply a GRanges object with an
#' additional column containing the copy number status (as integers) and
#' optionally another column containing the LOH status of each region.
#'
#' If there's only one label, it will be used to label both Copy Number and
#' LOH. If there are at least two labels, the first one is used to label the
#' Copy Number and the second one to labvel the LOH.
#'
#' If the function is called with a list of GRanges it plots every GRanges
#' independently as different tracks one below the other. Track positioning is
#' based on `autotrack` and the margin between track is controlled by
#' `track.margin`.
#'
#'
#' @usage plotCopyNumberCalls(karyoplot, cn.calls, cn.values=NULL, cn.column="cn", cn.colors=NULL,
#'                            loh.values=NULL, loh.column="loh", loh.color="#1E90FF", loh.height=0.3,
#'                            labels="", lab.cex=1, lab2.cex=NULL,
#'                            track.margin=0.01, r0=0, r1=1, ...)
#'
#' @param karyoplot A karyoplote object
#' @param cn.calls The CN calls to plot
#' @param cn.values The CN values. If NULL, they will be extracted from the cn.calls (defaults to NULL)
#' @param cn.column (defaults to "cn")
#' @param cn.colors (defaults to NULL)
#' @param loh.values The CN values. If NULL, they will be extracted from the cn.calls (defaults to NULL)
#' @param loh.column (defaults to "loh")
#' @param loh.color (defaults to "#1E90FF")
#' @param loh.height (defaults to 0.3)
#' @param labels (defaults to "")
#' @param lab.cex (defaults to 1)
#' @param lab2.cex (defaults to NULL)
#' @param track.margin (deafults to 0.01)
#' @param r0 (defaults to 0)
#' @param r1 (defaults to 1)
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





plotCopyNumberCalls <- function(karyoplot, cn.calls, cn.values=NULL, cn.column="cn", cn.colors=NULL,
                                loh.values=NULL, loh.column="loh", loh.color="#1E90FF", loh.height=0.3,
                                labels="", lab.cex=1, lab2.cex=NULL,
                                track.margin=0.01, r0=0, r1=1, ...) {

  #If cn.calls is a list, call this same function with each single element to actually produce the plot. Use autotrack to set the appropiate r0 and r1 values.
  if(methods::is(cn.calls, "list")) {
    for(i in seq_len(length(cn.calls))) {
      if(methods::is(labels, "list")) labels <- as.character(unlist(labels)) #If labels are a list, make them a vector
      if(length(labels)==length(cn.calls)) { #If there are as many labels as samples, assume each label should be used for one track
        lab <- labels[i]
      } else {
        lab <- labels[1]
      }
      rr <- karyoploteR::autotrack(current.track=i, total.tracks=length(cn.calls), margin = track.margin, r0=r0, r1=r1)
      plotCopyNumberCalls(karyoplot, cn.calls[[i]], r0=rr$r0, r1=rr$r1,
                          cn.values=cn.values, loh.values=loh.values,
                          cn.column=cn.column, loh.column = loh.column,
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


  if(is.null(cn.values)) {
    if(length(mcols(cn.calls))==0) stop("No cn.values given and cn.calls has no associated copy number data")
    #If no name for the copy number column was specified, use the first one
    if(is.null(cn.column)) {
      cn.values <- names(mcols(cn.calls))[1]
    } else {
      if(!(cn.column %in% names(mcols(cn.calls)))) stop("The cn.calls object does not have a column ", cn.column, ". No copy number data is available")
      cn.values <- mcols(cn.calls)[, cn.column]
    }
  }

  seg.cols <- segment.colors[as.character(cn.values)]


  karyoploteR::kpRect(karyoplot, data=cn.calls, y0=loh.height, y1=1, col=seg.cols, r0=r0, r1=r1, border=NA, ...)
  if("loh" %in% names(GenomicRanges::mcols(cn.calls))) {
    #convert loh=NA to no LOH
    cn.calls$loh[is.na(cn.calls$loh)] <- FALSE
    karyoploteR::kpRect(karyoplot, data=cn.calls[cn.calls$loh==TRUE], y0=0, y1=loh.height, r0=r0, r1=r1, col=loh.color, border=NA, ...)
  }
  if(!is.null(labels) && !is.na(labels) && all(is.character(labels) && length(labels)>0)) {
    if(length(labels)==1) {
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], r0=r0, r1=r1, cex=lab.cex, ...)
    } else { #If there are two labels, use the first one for the CN track an the second one for the LOH track
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], r0=r0+(r1-r0)*loh.height, r1=r1, cex=lab.cex, ...)
      karyoploteR::kpAddLabels(karyoplot, labels = labels[2], r0=r0, r1=r0+(r1-r0)*loh.height, cex=lab2.cex, ...)
    }
  }
  invisible(karyoplot)
}

