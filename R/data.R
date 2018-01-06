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

