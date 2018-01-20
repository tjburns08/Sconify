# Date: Jan 20, 2018
# Procedure: test the functions in the process.files script
# Purpose: unit testing for Sconify package

library(Sconify)
library(testthat)
context("Test the file processing arm of the Sconify package")


basal.file <- system.file('extdata',
                      'Bendall_et_al_Cell_Sample_C_basal.fcs',
                      package = "Sconify")

stim.file <- system.file('extdata',
                         'Bendall_et_al_Cell_Sample_C_IL7.fcs',
                         package = "Sconify")

markers.file <- system.file('extdata',
                        'markers.csv',
                        package = "Sconify")

input <- parse.markers(markers.file)[[1]]

# parse_markers
test_that("Markers csv file is successfully imported", {
    markers <- parse.markers(markers.file)
    expect_equal(length(markers), 2)
    expect_equal(length(markers[[1]]), 27)
    expect_equal(length(markers[[2]]), 16)
    expect_equal(names(markers), c("input", "functional"))
    expect_true(is.list(markers))
})

# fcs.to.tibble
test_that("Fcs file gets converted into a tibble data structure", {
    dat <- fcs.to.tibble(file = basal.file)
    expect_true(is.tibble(dat))
    expect_true(is.data.frame(dat))
    expect_false(is.matrix(dat))
    expect_true(is.atomic(dat[[1]]))
})

test_that("The asinh transform command works", {
    not.tr <- fcs.to.tibble(file = basal.file, transform = "none")
    tr <- fcs.to.tibble(file = basal.file)
    tr2 <- fcs.to.tibble(file = basal.file, transform = "asinh")
    expect_equal(tr[[3]], asinh(not.tr[[3]]/5))
    expect_equal(tr[[3]], tr2[[3]])
})

# Process.multiple.files
test_that("Processing multiple files works on a single file", {
    dat1 <- fcs.to.tibble(file = basal.file)
    dat2 <- process.multiple.files(files = basal.file, input = input)
    expect_equal(dat1[,input], dat2[,input])
    expect_equal(length(dat1) + 1, length(dat2))
    expect_equal(length(unique(dat2[["condition"]])), 1)
})

test_that("Process multiple files effectively sub-samples", {
    cell.number <- nrow(fcs.to.tibble(basal.file))
    testing <- c(cell.number,
                 cell.number %/% 2,
                 cell.number %/% 4,
                 cell.number %/% 8)

    lapply(testing, function(i) {
        curr <- process.multiple.files(files = basal.file, numcells = i,
                                       input = input)
        expect_equal(nrow(curr), i)
    })

    expect_error(process.multiple.files(files = basal.file, numcells = 0,
                                        input = input))

    expect_error(process.multiple.files(files = basal.file, numcells = -3,
                                        input = input))

    expect_error(process.multiple.files(files = basal.file, numcells = 10.76,
                                        input = input))

    expect_error(process.multiple.files(files = c(basal.file, stim.file), numcells = 1,
                                        input = input))
})

test_that("Process multiple files divdes the contribution of each file equally", {

    testing <- c(100, 99, 2)
    lapply(testing, function(i) {
        dat <- process.multiple.files(files = basal.file, numcells = 99,
                                      input = input)
        cond.dat <- dat[["condition"]]
        conds <- unique(cond.dat)
        cond1 <- cond.dat[cond.dat == conds[1]]
        cond2 <- cond.dat[cond.dat == conds[2]]
        expect_equal(length(cond1), length(cond2))
    })
})



