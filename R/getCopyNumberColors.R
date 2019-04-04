#' getCopyNumberColors
#'
#' @description
#' Get the colors representing the copy number levels
#'
#' @details
#' This function returns a vector of colors with names from 0 to the length
#' of the vector. The colors will be used to represent copy numbers from
#' 0 (homozygous deletion) to the copy number equal to the length of  the
#' vector (typically 7 or 8). The function has a single parameter `colors`,
#' which can be either the name of one of the available palettes or the
#' a vector of valid color specifications (i.e. c("red", "#FFAAAA", "gray",
#' "#AAAAFF", "#5555FF", "#0000FF")).  If `colors` is a palette name,
#' the palette will be returned. If color is a color vector, it will be returned
#' with the names set as needed.
#'
#' Available palettes are: "gree_orange_red" and "red_blue"
#'
#' @usage getCopyNumberColors(colors=NULL)
#'
#' @param colors  (character) Either a palette name (i.e. "green_orange_red") or a vector of colors. If NULL the default palette ("green_orange_red") will be returned. (defaults to NULL)
#'
#' @return
#' A vector of colors with the names set to the positions in the vector.
#'
#' @examples
#'
#' #getCopyNumberColors()
#'
#' #getCopyNumberColors("red_blue")
#'
#' #getCopyNumberColors(c("red", "#FFAAAA", "gray", "#AAAAFF", "#5555FF", "#0000FF"))
#'
#' @export getCopyNumberColors
#'
#'


getCopyNumberColors <- function(colors=NULL) {
  if(is.null(colors)) colors <- "green_orange_red"

  #Define the different palettes
  green_orange_red <- stats::setNames(c("darkgreen", "lightgreen", "gray","#FFD700", "#FFA500", "#EE7600", "#FF4500", "#CD0000", "#8B0000", "#0D0D0D"), as.character(0:8))
  red_blue <- stats::setNames(c("#EE0000", "#FFC1C1", "#E0E0E0", "#B2DFEE", "#87CEFA", "#1E90FF", "#0000FF"), c(0:6))

  if(is.character(colors) && length(colors)==1) {
    if(colors=="green_orange_red") {
      return(green_orange_red)
    } else {
      if(colors=="red_blue") {
        return(red_blue)
      }
    }
  } else {
    #assume colors is a valid
    if(all(karyoploteR::is.color(colors))) {
      return(stats::setNames(colors, seq_len(length(colors))-1))
    }
  }

  #If we arrive here, the given colors were not valid. Returning the default colors
  warning(paste0("Invalid colors. Returning default palette. Invalid: ", paste0(colors, collapse=", ")))
  return(green_orange_red)
}

