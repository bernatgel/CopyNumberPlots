library(CopyNumberPlots)
library(testthat)
context("Get Columns functions")


df <- data.frame("id"= "rs1234", "endogenous" = FALSE, "chromosome"="chr1", "Start"=0, "end.position"=100,
                 "copy.number.level"=3, "LOH"=0, "median.value.per.segment"=1.2,
                 "BAF"=0.2, "Log Ratio"=1.5,  "strange.name"="strange.value")

gr <- regioneR::toGRanges(data.frame(chr = "chr5", start = 102562, end = 224624, lrr = 0.5642))

#testing getColumn
testthat::test_that("getColumn returns an integer", {
  
  #test msg error. Either pattern or column must be provided
  expect_error(getColumn(df = df), "Either col or pattern must be provided")
  
  #getColumn returns an integer if pattern is given
  column <- getColumn(df = df, pattern = "", col = NULL, verbose = FALSE)
  expect_is(column, "integer")
  
  # if pattern is numeric
  
  #getColumn returns an integer if col is given
  column <- getColumn(df = df, pattern = NULL, col= 2, verbose = FALSE)
  expect_is(column, "integer")
})

testthat::test_that("col parameter is correct",{
    # length of col must be one
    expect_error(getColumn(df = df, pattern = NULL, col = c(2,1), verbose = FALSE),
                 "col parameter must be either NULL, a character of length 1 or a integer of length 1")
                    
    # col must have the same length as the length of df
    expect_error(getColumn(df = df, pattern = NULL, col = 55), "col must be a number between one and length of names of df")
    
    #true/false
    expect_error(getColumn(df= df, col = TRUE, verbose = FALSE),
                 "col parameter must be either NULL, a character of length 1 or a integer of length 1")
    #NULL
    expect_error(getColumn(df = df, col = NULL, verbose = FALSE),
                 "Either col or pattern must be provided")
    #NA
    expect_error(getColumn(df = df, col = NA, verbose = FALSE),
                 "col parameter must be either NULL, a character of length 1 or a integer of length 1")
    
    #granges
    expect_error(getColumn(df = df, col = gr, verbose = FALSE),
                 "col parameter must be either NULL, a character of length 1 or a integer of length 1")
    #col = -2
    expect_error(getColumn(df = df, col = -2, verbose = FALSE),
                 "col must be a single integer between one and length of names of df")
    #col = 2.5
    expect_error(getColumn(df = df, col = 2.5, verbose = FALSE),
                 "If col is numeric it must be an integer")
  
})

testthat::test_that("col.names is not NULL",{
  expect_error(getColumn(df = gr, col ="lrr"),
               "col.names cannot be NULL if we look for a pattern or column name.")
})

testthat::test_that(" testing msg.col.name parameter",{
  # msg.col.name must be a character
  expect_error(getColumn(df = df, pattern = "", msg.col.name = 125), "msg.col.name must be a character")
  expect_message(getColumn(df = df, pattern = "id", msg.col.name = NULL), "The column identified is: id")
  expect_message(getColumn(df = df, pattern = "id", msg.col.name = "ID"), "The column identified as ID is: id")

})


testthat::test_that("testing avoid pattern",{
  
  expect_equal(getColumn(df = df, pattern = "end", verbose = FALSE),2)
  expect_equal(getColumn(df = df, pattern = "end", avoid.pattern = NA, verbose = FALSE),2)
  expect_equal(getColumn(df = df, pattern = "end", avoid.pattern = NULL, verbose = FALSE),2)
  expect_equal(getColumn(df = df, pattern = "end", avoid.pattern = "endogenous", verbose = FALSE),5)

})


context("Testing transformChr")
testthat::test_that("testing if transformChr returs a character",{
  seg.file <- system.file("extdata", "DNACopy_output.seg", package = "CopyNumberPlots", mustWork = TRUE)
  seg.data <- read.table(file = seg.file, sep = "\t", skip = 1, stringsAsFactors = FALSE)
  colnames(seg.data) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  expect_is(transformChr(chr = seg.data$chrom, chr.transformation = "23:X"), "character")
  
  #testing chr.transformation equal to NULL
  expect_error(transformChr(chr = seg.data$chrom, chr.transformation = 1),
               "chr.transformation must be a character")
  expect_error(transformChr(chr = seg.data$chrom, chr.transformation = NULL),
               "chr.transformation parameter must be a character \"key:value\" of length one")
  expect_error(transformChr(chr = seg.data$chrom, chr.transformation = ""),
               "chr.transformation parameter must be a character \"key:value\" of length one")
})



context("Testing removeNAs")
testthat::test_that("testing removeNAs",{
  #testing if removeNAs return a GRanges or a list
  seg.data <- regioneR::toGRanges(data.frame(chr = c("chr1", "chr1", "chr2", "chr5"),
                                             start = c(0,50000,8014630,14523572),
                                             end = c(48953, 7023664,9216331,153245687),
                                             lrr = c(NA,0.25,1.5,NA)))
  expect_is(removeNAs(snp.data = seg.data, verbose = FALSE), "GRanges")
  
  seg.data <- list(a = regioneR::toGRanges(data.frame(chr = c("chr1", "chr1", "chr2", "chr5"), 
                                                           start = c(0,50000,8014630,14523572), 
                                                           end = c(48953, 7023664,9216331,153245687),
                                                           lrr = c(NA,0.25,1.5,NA),baf = c(1.5,2.5,NA,6),
                                                           id = c("rs52456","rs52457","rs52458","rs52459"))),
                        b=regioneR::toGRanges(data.frame(chr = c("chr1", "chr1", "chr2", "chr5"),
                                                         start = c(0,50000,8014630,14523572), 
                                                         end = c(48953, 7023664,9216331,153245687),
                                                         lrr = c(2.5,NA,1.5,0.25), baf = c(1.5,2.5,NA,6),
                                                         id = c("rs52456","rs52457","rs52458","rs52459"))))
  expect_is(removeNAs(snp.data = seg.data, verbose = FALSE), "list")
    
  #testing if we introduce a df instead of GRanges or a list of GRanges
  expect_error(removeNAs(snp.data = df), "snp.data must be a GRanges or a list of GRanges")
 
  #lrr.na is a character
  expect_error(removeNAs(snp.data = seg.data, lrr.na = "a"),
               "lrr.na, baf.na, id.na and verbose must be either TRUE or FALSE")
  
  
})

