#' Basal data
#'
#' The basal cells from a single patient in the Wanderlust dataset
#'
#' @format A tibble of 10,000 cells by 51 features. All markers in the
#' dataset, along with pre-calculated Wanderlust value and condition,
#' which is a string that denotes that this is the "basal" condition for
#' each row. Important when this is concatenated with additional conditions
"basal.data"

#' Bodenmiller-Zunder GM-CSF post-SCONE final data, that's been quantile
#' normalized and z scored
#'
#' The post-SCONE output from a per-marker quantile normalized and z scored
#  Bodenmiller-Zunder dataset pair of fcs files, one untreated and one
#' treated with GM-CSF.
#'
#' @format A tibble of 10,000 cells by 69 features. This includes all the
#' original parameters, the KNN-generated comparisons, differential
#' abundance ("fraction.cond.2), and two t-SNE coordinates.
"bz.gmcsf.final.norm.scale"

#' Bodenmiller-Zunder GM-CSF post-SCONE final data
#'
#' The post-SCONE output from Bodenmiller-Zunder dataset pair of fcs files,
#' one untreated and one treated with GM-CSF.
#'
#' @format A tibble of 9,998 cells by 69 features. This includes all the
#' original parameters, the KNN-generated comparisons, differential
#' abundance ("fraction.cond.2), and two t-SNE coordinates.
"bz.gmcsf.final"

#' The Wanderlust dataset, combined basal and IL7 conditions
#'
#' A single patient pair of basal and IL7 treated cells from bone marrow
#' gated for B cell precursors.
#'
#' @format A tibble of 10,000 cells by 51 features, including all the input
#' markers, Wanderlust values, and the condition. The first 5000 rows
#' are untreated cells and the last 5000 rows are IL7 treated.
"combined"

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
#' @format A tibble of 10,000 cells and 87 feaures, including the input
#' features, the SCONE-generated comparisons, differential abundance, and
#' two t-SNE dimesnions
"final"

#' A 2000 cell subset of the "final" data
#'
#' 2000 cells randomly taken from the 10,000 cell final dataset.
#'
#' @format A tibble of 2000 cells and 87 feaures, including the input
#' features, the SCONE-generated comparisons, differential abundance, and two
#' t-SNE runs. One genrated from the 10,000 cell data (bh-SNE1 and bh-SNE2),
#' along with a run from the 2000 cell data (bh-SNE11 and bh-SNE21). This
#' is used in the post-processing vignette to describe the subsample.and.tsne
#' function from this package
"final.sub"

#' Functional markers from the Wanderlust dataset.
#'
#' These are the markers that will be used in the KNN comparisons, as opposed
#' to the KNN generation.
#'
#' @format A vector of strings.
"funct.markers"

#' A named vector to help the user determine the ideal k.
#'
#' This is the output of the impute.testing function used on the
#' Wanderlust dataset, which finds the avergae imputation error of all signal
#' markers imputed from KNN of surface markers.
#'
#' @format A named vector, where the elements are averge imputation error
#' and the names are the values of from a 10,000 cell dataset.
"ideal.k"

#' IL7 data
#'
#' The IL7 treated cells from a single patient in the Wanderlust dataset
#'
#' @format A tibble of 10,000 cells by 51 features. All markers in the
#' dataset, along with pre-calculated Wanderlust value and condition,
#' which is a string that denotes that this is the "IL7" condition for
#' each row. Important when this is concatenated with additional conditions
"il7.data"

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

#' K titration
#'
#' A titration of values of K used for a 10,000 cell dataset as input
#' for the impute.testing function.
#'
#' @format A vector of numbers
"k.titration"

#' KNN density
#'
#' A vector of per-cell KNN density from the Wanderlust dataset, generated
#' from the get.knn.de command. Of note, the SCONE pipeline will automatically
#' generate this for you. This is just one divided by the average distance
#' to KNN.
#'
#' @format A numeric vector
"knn.density"

#' KNN list
#'
#' A list of the features of each cell's KNN. This is from the optional
#' make.knn.list function, that allows the user to do KNN-based operations
#' that go beyond the scope of the Sconify package.
#'
#' @format A list of tibbles, with each element being a tibble of cells
#' by features corresponding to a given cell's KNN. The element of the list
#' corresponds to the cell identity (row) of the original input data.
"knn.list"

#' Markers for the Wanderlust dataset
#'
#' Both the surface and functional markers for the Wanderlust dataset
#'
#' @format a tibble with two columns, "surface" and "fucntional."
"markers"

#' Multiple donor dataset post-SCONE
#'
#' The Fragidakis dataset with multiple donors, post-SCONE analysis.
#'
#' @format A tibble of 18,514 cells by 124 features. These features include
#' original input markers, per-donor comparisons, per-cell comparisons,
#' density, density estimation, and per-condition differential abundance.
"md.final"

#' Multiple donor dataset input markers
#'
#' The markers to be used for KNN generation for the Fragidakis multiple donor
#' dataset. These are primarily surface markers
#'
#' @format A vector of strings
"md.input"

#' Fragidakis multiple donor dataset nearest neighbor information
#'
#' Two matrices corresponding to the KNN identity and distance.
#'
#' @format List of 2. First element named nn.index is a matrix of 19,992 cells
#' by 200 nearest neighbor positions, with the matrix elements being cell
#' identity.The second element is a matrix of 19,992 cells by 200 nearest
#' neighborhood positions with the matrix elements being distances.
"md.nn"














