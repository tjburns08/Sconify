# Date: Janurary 26, 2018
# Procedure: Unit testing for the knn and statistics script
# Purpose: Unit testing for the Sconify package
library(testthat)
library(Sconify)
context("Test the knn and statistics arm of the sconify package")

k <- 50
test.nn <- fnn(wand.combined, input.markers = input.markers, k = k)
test.scone <- scone.values(nn.matrix = test.nn,
                    cell.data = wand.combined,
                    scone.markers = funct.markers,
                    unstim = "basal")

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

test_that("knn density estimation is produced", {
    tmp <- get.knn.de(test.nn)
    expect_equal(length(tmp), nrow(wand.combined))
    expect_true(all(tmp > 0))
})

test_that("knn density estimation requires the list of two from fnn output", {
    expect_error(get.knn.de(test.nn[[1]]))
    expect_error(get.knn.de(test.nn[[2]]))
})

test_that("knn list is created for each cell", {
    tmp <- make.knn.list(wand.combined, test.nn)
    expect_equal(length(tmp), nrow(wand.combined))
    expect_true(all(colnames(tmp[[1]]) == colnames(wand.combined)))
    expect_equal(nrow(tmp[[1]]), k)
})

test_that("Scone values outputs a tibble of statistical values", {
    expect_equal(ncol(test.scone), 2*length(funct.markers) + 2)
})

test_that("Scone values produces a proper knn density estimation", {
    expect_true(all(test.scone$density, get.knn.de(test.nn)))
})

test_that("Scone produces a proper readout of differential abundance", {
    expect_false(all(test.scone$IL7.fraction.cond.2 > 1))
})

test_that("Scone wont perform statistics unless a proper test name is used", {
    expect_error(scone.values(nn.matrix = test.nn,
                              scone.markers = funct.markers,
                              unstim = "basal",
                              stat.test = "tyler's test"))
})

test_that("Scone wont perform statistics unless a proper basal name is used", {
    expect_error(scone.values(nn.matrix = test.nn,
                              scone.markers = funct.markers,
                              unstim = "bas"))
})

test_that("Scone produces FDR adjusted q-values", {
    expect_true(all(test.scone$`Ki67(Sm152)Di.IL7.qvalue` <= 1))
    expect_false(all(test.scone$`Ki67(Sm152)Di.IL7.qvalue` < 1)) # p.adjust
})

test_that("Scone produces fold q thresholded fold changes", {
    tmp1 <- test.scone$`pSTAT5(Nd150)Di.IL7.qvalue` # 1.0
    tmp2 <- test.scone$`pSTAT5(Nd150)Di.IL7.change`
    expect_equal(which(tmp1 < 0.05), which(tmp2 > 0))
})








