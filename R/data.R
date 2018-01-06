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
#' The post-SCONE output from a normalized and z scored Bodenmiller-Zunder
#' dataset pair of fcs files, one untreated and one treated with GM-CSF.
#'
#' @format A tibble of 10,000 cells by 69 features. This includes all the
#' original parameters, the KNN-generated comparisons, differential
#' abundance ("fraction.cond.2), and two t-SNE coordinates.
"bz.gmcsf.final.norm.scale"
