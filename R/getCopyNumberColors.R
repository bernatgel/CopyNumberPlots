#' getCopyNumberColors
#'
#' @description
#' Get the colors representing each copy number
#'
#' @details
#' 
#'
#' @usage getCopyNumberColors()
#'
#' @param colors
#'
#'
#' @return
#' A vector with colors
#'
#' @examples
#'
#'
#'
#' @export getCopyNumberColors
#'
#'


getCopyNumberColors <- function(colors) {
  if(is.null(colors)) colors <- "green_orange_red"

  #Define the different palettes
  green_orange_red <- setNames(c("darkgreen", "lightgreen", "gray","#FFD700", "#FFA500", "#EE7600", "#FF4500", "#CD0000", "#8B0000", "#0D0D0D"), as.character(0:8))
  
  if(is.character(colors) && length(colors)==1) {
    if(colors=="green_orange_red") {
      return(green_orange_red)
    }
  } else {
    #assume colors is a valid
    if(all(is.color(colors))) {
      return(setNames(colors, seq_len(length(colors))))
    }
  }
  
  #If we arrive here, the given colors were not valid. Returning the default colors
  warning(paste0("Invalid colors. Returning default palette. Invalid: ", paste0(colors, collapse=", ")))
  return(green_orange_red)
}

#From the "network" package
is.color<-function(x){
  xic<-rep(FALSE,length(x))         #Assume not a color by default
  xc<-sapply(x,is.character)        #Must be a character string
  #For characters, must be a named color or a #RRGGBB/#RRGGBBAA sequence
  xic[xc]<-(x[xc]%in%colors()) | ((nchar(x[xc])%in%c(7,9))&(substr(x[xc],1,1)=="#"))
  xic[is.na(x)]<-NA                 #Missing counts as missing
  #Return the result
  xic
}
