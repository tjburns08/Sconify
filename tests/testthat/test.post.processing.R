# Date: January 28, 2018
# Procedure: testing the post processing function
# Purpose: unit testing for the testthat package
library(testthat)
library(Sconify)
context("Test the post processing functions for the Sconify package")

test_that("Subsampling and tSNE functionality works", {
    test.tsne <- subsample.and.tsne(wand.combined, input.markers, 100)
    expect_equal(ncol(test.tsne), ncol(wand.combined) + 2)
    expect_equal(length(grep("bh-SNE", colnames(test.tsne))), 2)
})

test_that("String to numbers works", {
    test.str <- string.to.numbers(c("hi", "hi", "bye", "hi"))
    expect_equal(test.str, c(1, 1, 2, 1))
})

test.pp <- post.processing(scone.output = wand.scone, cell.data = wand.combined, input = input.markers, tsne = TRUE, log.transform.qvalue = TRUE)
test_that("Basic post-processing functionality", {
    expect_equal(ncol(test.pp), ncol(wand.final))
})
