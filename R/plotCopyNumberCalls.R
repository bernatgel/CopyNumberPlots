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
#' @usage plotCopyNumberCalls(karyoplot, cn.calls, cn.values=NULL, cn.column="cn", cn.colors=NULL, loh.values=NULL, loh.column="loh", loh.color="#1E90FF", loh.height=0.3, labels=NULL, label.cex=1, label2.cex=NULL,track.margin=0.01, r0=0, r1=1, ...)
#'
#' @param karyoplot A karyoplote object
#' @param cn.calls (a GRanges, a list of GRanges or a GRangesList) An object with the positions of the CN calls and a column with the CN values. Other columns are ignored. If it's a list of GRanges with different samples, all samples will be plotted, splitting the total plot space between them.
#' @param cn.values (integer vector) The CN values. If NULL, they will be extracted from the cn.calls (defaults to NULL)
#' @param cn.column (integer or character vector) The name or number of the column with CN information.(defaults to "cn")
#' @param cn.colors (colors) The colors assigned to gains and losses (defaults to NULL)
#' @param loh.values (logical vector) A logical vector that indicates whether there is or not LOH or a vector that can be coerced as logical. (defaults to NULL)
#' @param loh.column (number or character) The name or number of the column with LRR information. (defaults to "loh")
#' @param loh.color (a color) The color assigned to LOH values. (defaults to "#1E90FF")
#' @param loh.height (numeric) The proportion of r0 and r1 of the vertical space over each chromosome dedicated to loh. It is dedicated the 30\% of the vertical space by default.(defaults to 0.3)
#' @param labels (character) The text of the label to identify the data. If NA, no label will be plotted. If NULL, if snps is a single sample GRanges it will default to "CN", if it's a list of samples it will default to the names in the list or consecutive numbers if names(snps) is NULL. (defaults to NULL)
#' @param label.cex (numeric) The size of tthe label (defaults to 1)
#' @param label2.cex (numeric) The size of the label 2. If NULL label2.cex will be lable.cex. (defaults to NULL)
#' @param track.margin (numeric) If cn.calls is a list object, this is the margin between the samples CN. (deafults to 0.01)
#' @param r0 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 0)
#' @param r1 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 1)
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
#' s1.calls.file <- system.file("extdata", "S1.segments.txt", package = "CopyNumberPlots", mustWork = TRUE)
#' s1.calls <- loadCopyNumberCalls(s1.calls.file)
#'
#' s2.calls.file <- system.file("extdata", "S2.segments.txt", package = "CopyNumberPlots", mustWork = TRUE)
#' s2.calls <- loadCopyNumberCalls(s2.calls.file)
#'
#' kp <- plotKaryotype(chromosomes="chr1")
#' plotCopyNumberCalls(kp, s1.calls)
#'
#' kp <- plotKaryotype(chromosomes="chr1")
#' plotCopyNumberCalls(kp, s1.calls, cn.colors="red_blue")
#' 
#' #List of GRanges
#' cn.calls <- list(s1=s1.calls, s2 =s2.calls)
#' 
#' kp <-plotKaryotype(chromosomes="chr1")
#' plotCopyNumberCalls(kp, cn.calls, cn.colors="red_blue")
#'
#' @export plotCopyNumberCalls
#'
#' @import karyoploteR
#' @importFrom IRanges subsetByOverlaps
#'
plotCopyNumberCalls <- function(karyoplot, cn.calls, cn.values=NULL, cn.column="cn", cn.colors=NULL,
                                loh.values=NULL, loh.column="loh", loh.color="#1E90FF", loh.height=0.3,
                                labels=NULL, label.cex=1, label2.cex=NULL,
                                track.margin=0.01, r0=0, r1=1, ...) {

  #If cn.calls is a list, call this same function with each single element to actually produce the plot. Use autotrack to set the appropiate r0 and r1 values.
  if(methods::is(cn.calls, "list") || methods::is(cn.calls, "CompressedGRangesList")) {
    labels <- prepareLabels(labels = labels, x = cn.calls)
    for(i in seq_len(length(cn.calls))) {
      rr <- karyoploteR::autotrack(current.track=i, total.tracks=length(cn.calls), margin = track.margin, r0=r0, r1=r1)
      plotCopyNumberCalls(karyoplot, cn.calls[[i]], r0=rr$r0, r1=rr$r1,
                          cn.values=cn.values, loh.values=loh.values,
                          cn.column=cn.column, loh.column = loh.column,
                          cn.colors = cn.colors, loh.color = loh.color,
                          labels = labels[i], label.cex = label.cex, label2.cex = label2.cex,
                          loh.height = loh.height, ...)
    }
    return(invisible(karyoplot))
  }
  
  
  if(is.null(label2.cex)) label2.cex <- label.cex
  
  #use toGRanges to build a GRanges if cn.calls was anything else
  cn.calls <- regioneR::toGRanges(cn.calls)
  
  if(is.null(labels)) labels <- "CN"
  
  if(!is.null(labels) && !is.na(labels) && all(is.character(labels) && length(labels)>0)) {
    if(length(labels)==1) {
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], r0=r0, r1=r1, cex=label.cex, ...)
    } else { #If there are two labels, use the first one for the CN track an the second one for the LOH track
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], r0=r0+(r1-r0)*loh.height, r1=r1, cex=label.cex, ...)
      karyoploteR::kpAddLabels(karyoplot, labels = labels[2], r0=r0, r1=r0+(r1-r0)*loh.height, cex=label2.cex, ...)
    }
  }
  
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
  
  # loh values
  if(is.null(loh.values)) {
    if(loh.column %in% names(GenomicRanges::mcols(cn.calls))){
      loh.values <- cn.calls[,loh.column]
     }
  }
  
  if (!is.null(loh.values)){
    if(!is.logical(loh.values)) loh.values <- tryCatch(as.logical(loh.values))
    #convert loh=NA to no LOH
    loh.values[is.na(loh.values)] <- FALSE
  }
  

  
  segment.colors <- getCopyNumberColors(cn.colors)

  #Explicitly filter the segments, since it will use the wrong colors otherwise
  #Check: is it still needed? The bug should be fixed
  plot.region <- karyoplot$plot.region
  segments <- IRanges::subsetByOverlaps(cn.calls, plot.region)

  seg.cols <- segment.colors[as.character(cn.values)]

  karyoploteR::kpRect(karyoplot, data=cn.calls, y0=loh.height, y1=1, col=seg.cols, r0=r0, r1=r1, border=NA, ...)
  
  if(!is.null(loh.values)){
    karyoploteR::kpRect(karyoplot, data=cn.calls[loh.values], y0=0, y1=loh.height, r0=r0, r1=r1, col=loh.color, border=NA, ...)
  }
 
  
  invisible(karyoplot)
}

