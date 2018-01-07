#' Bodenmiller-Zunder GM-CSF post-SCONE final data, that's been quantile
#' normalized and z scored.
#'
#' The post-SCONE output from a per-marker quantile normalized and z scored
#  Bodenmiller-Zunder dataset pair of fcs files, one untreated and one
#' treated with GM-CSF. We ran this on 10,000 cells and sub-sampled
#' to 1000 for this package.
#'
#' @format A tibble of 1000 cells by 69 features. This includes all the
#' original parameters, the KNN-generated comparisons, differential
#' abundance ("fraction.cond.2), and two t-SNE coordinates.
"bz.gmcsf.final.norm.scale"

#' Bodenmiller-Zunder GM-CSF post-SCONE final data
#'
#' The post-SCONE output from Bodenmiller-Zunder dataset pair of fcs files,
#' one untreated and one treated with GM-CSF.We ran this on 10,000 cells and
#' subsampled to 1000 for this vignette.
#'
#' @format A tibble of 1000 cells by 69 features. This includes all the
#' original parameters, the KNN-generated comparisons, differential
#' abundance ("fraction.cond.2), and two t-SNE coordinates.
"bz.gmcsf.final"

#' Wanderlust data combined basal and IL7 cells
#'
#' A single patient pair of basal and IL7 treated cells from bone marrow
#' gated for B cell precursors.
#'
#' @format A tibble of 1000 cells by 51 features, including all the input
#' markers, Wanderlust values, and the condition. The first 500 rows
#' are untreated cells and the last 500 rows are IL7 treated.
"wand.combined"

#' Random musing
#'
#' Seriously random
#'
#' @format a string
"exist"

#' Post-scone output of the "combiend" Wanderlust data.
#'
#' "combined" data taken through KNN generation and comparisons, along
#' with t-SNE map generation.
#'
#' @format A tibble of 1000 cells and 87 feaures, including the input
#' features, the SCONE-generated comparisons, differential abundance, and
#' two t-SNE dimesnions
"wand.final"

#' Functional markers from the Wanderlust dataset.
#'
#' These are the markers that will be used in the KNN comparisons, as opposed
#' to the KNN generation.
#'
#' @format A vector of strings.
"funct.markers"

#' A named vector to help the user determine the ideal k
#' for the Wanderlust dataset.
#'
#' This is the output of the impute.testing function used on the
#' Wanderlust dataset, which finds the avergae imputation error of all signal
#' markers imputed from KNN of surface markers.
#'
#' @format A named vector, where the elements are averge imputation error
#' and the names are the values of from a 10,000 cell dataset.
"wand.ideal.k"

#' Wanderlust IL7 data
#'
#' The IL7 treated cells from a single patient in the Wanderlust dataset
#'
#' @format A tibble of 1000 cells by 51 features. All markers in the
#' dataset, along with pre-calculated Wanderlust value and condition,
#' which is a string that denotes that this is the "IL7" condition for
#' each row. Important when this is concatenated with additional conditions
"wand.il7"

#' Input markers for the Wanderlust dataset
#'
#' These are the markers that KNN generation will be done on for the
#' Wanderlust dataset. These are mostly surface markers.These are the
#' same markers one would use as input for clustering or t-SNE generation,
#' for exmaple, as they are not expected to change through the duration
#' of the quick IL7 stimulation.
#'
#' @format A vector of strings corresponding to the markers.
"input.markers"

#' Markers for the Wanderlust dataset
#'
#' Both the surface and functional markers for the Wanderlust dataset
#'
#' @format a tibble with two columns, "surface" and "fucntional."
"markers"

#' Wanderlust scone output
#'
#' The scone output for the Wanderlust dataset
#'
#' @format A tibble of 1000 cells by 34 features. These features include
#' the KNN comparisons, KNN density estimation, and differential abundance.
#' Note that this tibble gets concatenated with the original tibble, as well
#' as two t-SNE dimensions in the post.processing() command of the pipeline.
"wand.scone"










