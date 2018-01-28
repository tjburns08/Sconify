# Date: January 28, 2018
# Procedure: testing the k selection
# Purpose: unit testing
context("Testing the functionality of the ideal k finder")

k.titration <- c(5, 50, 500)
test <- impute.testing(k.titration = k.titration,
                       cells = wand.il7,
                       input.markers = input.markers,
                       test.markers = funct.markers)

test_that("High level desired output from the impute testing function", {
    expect_equal(length(test), length(k.titration))
})

test_that("Convex loss function expected", {
    expect_true(test[1] > test[2])
    expect_true(test[2] < test[3])
})

test_that("Edge case", {
    expect_error()
})
