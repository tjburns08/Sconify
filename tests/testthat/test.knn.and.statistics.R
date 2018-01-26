# Date: Janurary 26, 2018
# Procedure: Unit testing for the knn and statistics script
# Purpose: Unit testing for the Sconify package
library(testthat)
library(Sconify)
context("Test the knn and statistics arm of the sconify package")

k <- 100
test.nn <- fnn(wand.combined, input.markers = input.markers, k = k)

test_that("fnn function produces a list of two", {
    expect_equal(length(test.nn), 2)
})

test_that("fnn produces a matrix of n rows (number of cells) and k columns (
          number of nearest neighbors", {
    tmp <- test.nn[[1]]
    expect_equal(ncol(tmp), k)
    expect_equal(nrow(tmp), nrow(wand.combined))
})

test_that("fnn does not work when there is no input markers", {
    expect_error(fnn(wand.combined, input.markers = c("hi"), k = k))
    expect_error(fnn(wand.combined, input.markers = c(), k = k))
    expect_error(fnn(wand.combined, input.markers = c("CD3"), k = k))
})

test_that("fnn does not work with some k values", {
    expect_error(fnn(wand.combined, input.markers = input.markers, k = 0))
    expect_error(fnn(wand.combined, input.markers = input.markers, k = -3))
})

test_that("fnn ")
