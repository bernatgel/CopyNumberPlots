# Data frame parsing functions

#' getColumn
#'
#' @description
#' Use simple pattern matching to try to identify a column in a data.frame by 
#' its name
#'
#' @details
#' This function will use pattern matching to try to identify which column of 
#' a data.frame contains a certain data. It will return its number. If 
#' \code{col} is specified it will not use pattern matching but identify
#' if a column with that exact name exists and return its position.
#'
#'
#' @usage getColumn(df, col=NULL, pattern="",  msg.col.name="", needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, teh column will be searched using pattern (defaults to NULL)
#' @param pattern (character) The pattern to match the column name. If more than one column matches the pattern, the first one (leftmost) will be returned. The pattern may be any valid regular expression. (defaults to "")
#' @param msg.col.name (character) Only used in the error message to make the message clearer. The name of the column we are searching for. (defaults to "")
#' @param needed (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#'
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#'
#' getColumn(df, pattern="Chr|chr",  msg.col.name="Chromosome", needed=TRUE)
#' getColumn(df, col="chromosome",  msg.col.name="Chromosome", needed=TRUE)
#'
#' @export getColumn



getColumn <- function(df, col=NULL, pattern="",  msg.col.name="", needed=TRUE) {
  col.num <- integer(0)
  if(is.null(col)) {
    if(length(pattern)>0) {
      col.num <- which(grepl(names(df), pattern = pattern, ignore.case = TRUE))[1]
    } else {
      stop("Either col or pattern must be provided")
    }
  } else {
    if(is.numeric(col)) {
      col.num <- col
    } else {
      if(is.character(col))  {
        col.num <- which(names(df)==col)
      }
    }
  }
  if(is.na(col.num) || length(col.num)==0) {
    if(needed==TRUE) {
      stop("The column ", msg.col.name, " was not found in the data")
    } else {
      col.num <- NULL
    }
  }
  return(col.num)
}



#' getChrColumn
#'
#' @description
#' Identify the column in a data.frame with the chromosome information
#'
#' @details
#' Identify the column of a data.frame that contains the chromosome 
#' information and return its position
#'
#' @usage getChrColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getChrColumn(df)
#' getChrColumn(df, "strange.name")
#'
#' @export getChrColumn
getChrColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="chr", msg.col.name="Chromosome", needed=needed))
}

#' getPosColumn
#'
#' @description
#' Identify the column in a data.frame with the position information
#'
#' @details
#' Identify the column of a data.frame that contains the position 
#' information and return its position
#'
#' @usage getPosColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#'
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getPosColumn(df)
#' getPosColumn(df, "strange.name")
#'
#' @export getPosColumn
getPosColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="Position|pos|loc|maploc", msg.col.name = "Position", needed=needed))
}

#' getStartColumn
#'
#' @description
#' Identify the column in a data.frame with the start position information
#'
#' @details
#' Identify the column of a data.frame that contains the position 
#' information and return its position
#'
#' @usage getStartColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getStartColumn(df)
#' getStartColumn(df, "strange.name")
#'
#' @export getStartColumn
getStartColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="start|first|begin", msg.col.name = "Start", needed=needed))
}

#' getEndColumn
#'
#' @description
#' Identify the column in a data.frame with the end position information
#'
#' @details
#' Identify the column of a data.frame that contains the position 
#' information and return its position
#'
#' @usage getEndColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#'
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getEndColumn(df)
#' getEndColumn(df, "strange.name")
#'
#' @export getEndColumn
getEndColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="end|last", msg.col.name = "End", needed=needed))
}


#' getCopyNumberColumn
#'
#' @description
#' Identify the column in a data.frame with the copy number information
#'
#' @details
#' Identify the column of a data.frame that contains the copy number
#' information and return its position
#'
#' @usage getCopyNumberColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#'
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getCopyNumberColumn(df)
#' getCopyNumberColumn(df, "strange.name")
#'
#' @export getCopyNumberColumn
getCopyNumberColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="cn|copy", msg.col.name = "Copy Number", needed=needed))
}

#' getLOHColumn
#'
#' @description
#' Identify the column in a data.frame with LOH information
#'
#' @details
#' Identify the column of a data.frame that contains the LOH 
#' information and return its position
#'
#' @usage getLOHColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#'
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getLOHColumn(df)
#' getLOHColumn(df, "strange.name")
#'
#' @export getLOHColumn
getLOHColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="loh|loss", msg.col.name = "LOH", needed=needed))
}

#' getSegmentValueColumn
#'
#' @description
#' Identify the column in a data.frame with the position of the segment value information
#'
#' @details
#' Identify the column of a data.frame that contains the segment 
#' information and return its position
#'
#' @usage getSegmentValueColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' 
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getSegmentValueColumn(df)
#' getSegmentValueColumn(df, "strange.name")
#'
#' @export getSegmentValueColumn
getSegmentValueColumn <- function(df, col=NULL, needed=TRUE) {
  return(getColumn(df, col=col, pattern="value|mean|median|ratio", msg.col.name = "Segment Value", needed=needed))
}

#' getBAFColumn
#'
#' @description
#' Identify the column in a data.frame with the bi-allelic frecuency information
#'
#' @details
#' Identify the column of a data.frame that contains the bi-allelic frecuency 
#' information and return its position
#'
#' @usage getBAFColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' 
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getBAFColumn(df)
#' getBAFColumn(df, "strange.name")
#'
#' @export getBAFColumn

getBAFColumn <- function(df, col = NULL, needed = TRUE){
  return(getColumn(df, col = col, pattern = "BAF|B.Allele|freq", msg.col.name = "B-Allele Frequency", needed = needed))
} 

#' getLRRColumn
#'
#' @description
#' Identify the column in a data.frame with the Log Ratio information
#'
#' @details
#' Identify the column of a data.frame that contains the Log Ratio 
#' information and return its position
#'
#' @usage getLRRColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' 
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getLRRColumn(df)
#' getLRRColumn(df, "strange.name")
#'
#' @export getLRRColumn

getLRRColumn <- function(df, col = NULL, needed = TRUE){
  return(getColumn(df, col = col, pattern = "LRR|Log.R.Ratio|Log", msg.col.name = "Log Ratio", needed = needed))
} 

#' getIDColumn
#'
#' @description
#' Identify the column in a data.frame with the ID information
#'
#' @details
#' Identify the column of a data.frame that contains ID 
#' information and return its position
#'
#' @usage getIDColumn(df, col=NULL, needed=TRUE)
#'
#' @param df (data.frame or equivalent) The data frame were columns are searched
#' @param col (number of character) If a number, it will be returned. If a character, it will be treated as a the exact name of the column. If a column with that name exists, its position in the data.frame will be returned. If it does not, NULL will be returned and an error might be raised (depending on the value of \code{needed}). If NULL, the column will be searched using a predefined pattern (defaults to NULL)
#' @param needed Whether the column is needed or not. If TRUE, an error will be raised if the column is not found. (defaults to TRUE)
#' 
#' @return
#' The number of the column matching the specification or NULL if no column was found.
#'
#' @examples
#' 
#' df <- data.frame("id"= "rs1234","chromosome"="chr1", "Start"=0, "end.position"=100,
#' "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
#' "BAF"=0.2, "Log Ratio"=1.5, "strange.name"="strange.value")
#' getIDColumn(df)
#' getIDColumn(df, "strange.name")
#'
#' @export getIDColumn

getIDColumn <- function(df, col = NULL, needed = TRUE){
  return(getColumn(df, col = col, pattern = "name|id|snp", msg.col.name = "Identifier", needed = needed))
} 


