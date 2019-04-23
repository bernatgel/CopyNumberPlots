library(CopyNumberPlots)
context("Get Columns functions")

#testing getColumn
testthat::test_that("getColumn returns an integer", {
  df <- data.frame("id"= "rs1234", "endogenous" = FALSE, "chromosome"="chr1", "Start"=0, "end.position"=100,
              "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
              "BAF"=0.2, "Log Ratio"=1.5,  "strange.name"="strange.value")
  
  #test msg error. Either pattern or column must be provided
  expect_error(getColumn(df = df), "Either col or pattern must be provided")
  
  #getColumn returns an integer if pattern is given
  column <- getColumn(df = df, pattern = "", col = NULL)
  expect_is(column, "integer")
  
  # if pattern is numeric
  
  #getColumn returns an integer if col is given
  column <- getColumn(df = df, pattern = NULL, col= 2)
  expect_is(column, "integer")
  
  # length of col must be one
  expect_error(getColumn(df = df, pattern = NULL, col = c(2,1)),
               "col parameter must be either NULL, a character of length 1 or a integer of length 1")

  
})


col.num <- getColumn(df = df, pattern = "Chr|chr",  msg.col.name = "Chromosome", needed = TRUE)
class(col.num)                     
