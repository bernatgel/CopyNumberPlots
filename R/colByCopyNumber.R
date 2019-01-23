#' colByCopyNumber
#'
#' @description
#' Assign a color to each data element according to the copy number status
#' of the its genomic location
#'
#' @details
#' Given a set of data elements positioned on the genome (either a GRanges
#' or any other format accepted by regioneR's \code{\link[regioneR]{toGRanges}})
#' and a set of copy number calls, assign a color to each element depending on
#' the copy number status of the genomic region where it's located. The colors
#' can be customized using the \code{cn.colors} parameter as documented in
#' \link{getCopyNumberColors}.
#'
#' @usage colByCopyNumber(data, cn.calls, cn.column="cn", cn.colors=NULL)
#'
#' @param data The data elements. Can be a GRanges or anything accepted by \code{\link[regioneR]{toGRanges}}
#' @param cn.calls The CN calls to use to select the color for the data elements
#' @param cn.column The name of the column with the copy number information. (defaults to "cn")
#' @param cn.colors The specification of the colors. Internally, it used to call the function \link{getCopyNumberColors}. Detailed documentation on color specification can be found there. (defaults to NULL)
#'
#'
#' @return
#' Return a vector of colors with the same length as the number of elements in
#' data
#'
#' @examples
#'
#' pos <- sort(floor(runif(1000, 1, 10000)))
#' baf.data <- toGRanges("chr1", pos, pos)
#' baf.data$baf <- rnorm(1000, mean = 0.5, sd = 0.05)
#' baf.data$baf[1:400] <- baf.data$baf[1:400] + c(0.2, -0.2)
#'
#' breakpoint <- start(baf.data)[400]
#' cn.calls <- toGRanges(data.frame(chr=c("chr1", "chr1"), start=c(1, breakpoint+1), end=c(breakpoint, 10000), cn=c(3,2)))
#'
#' kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
#' plotBAF(kp, baf.data, points.col = colByCopyNumber(baf.data, cn.calls))
#'
#'
#' @seealso \link[regioneR]{toGRanges}, \link{getCopyNumberColors}, \link{plotBAF}
#'
#' @export colByCopyNumber
#'
#'
#'


colByCopyNumber <- function(data, cn.calls, cn.column="cn", cn.colors=NULL) {

  segment.colors <- getCopyNumberColors(cn.colors)

  reg.cols <- as.character(segment.colors[as.character(as.data.frame(mcols(cn.calls)[cn.column])[,1])])

  cols <- karyoploteR::colByRegion(data, regions = cn.calls, colors = reg.cols)

  return(cols)
}

