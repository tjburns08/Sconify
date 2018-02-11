#' @import rflann
#' @importFrom stats dist
NULL

#' @title Imputes values for all markers (used as input) for each cell
#'
#' @description This function takes as input the markers to be imputed from
#' a pre-existing KNN computation.
#'
#' @param cells the input matrix of cells
#' @param input.markers the markers the user wants to impute
#' @param nn the matrix of k-nearest neighbors
#' (derived perhaps NOT from the "input markers" above)
#' @return a data frame of imputed cells for the "input markers" of interest
Impute <- function(cells, input.markers, nn) {
    result <- lapply(seq_len(nrow(cells)), function(i) {
        curr.nn <- cells[nn[i,],][,input.markers] %>% apply(., 2, mean)
    })

    result <- do.call(rbind, result) %>% as.tibble()
    return(result)
}

#' @title Impute testing
#'
#' @description Tests the euclidean distance error for imputation using knn
#' and markers of interest
#'
#' @param k.titration a vector integer values of k to be tested
#' @param cells a matrix of cells by features used as original input
#' @param input.markers markers to be used for the knn calculation
#' @param test.markers the markers to be tested for imputation
#' (either surface or scone)
#' @return the median imputation error for each value k tested
#' @examples
#' ImputeTesting(k.titration = c(10, 20),
#'   cells = wand.combined,
#'   input.markers = input.markers,
#'   test.markers = funct.markers)
#' @export
ImputeTesting <- function(k.titration, cells, input.markers, test.markers) {
    final.distances <- lapply(k.titration, function(k) {
        nn <- Fnn(cell.df = cells, input.markers = input.markers, k = k)[[1]]
        imputed.cells <- Impute(cells = cells, input.markers = test.markers,
                                nn = nn)

        print(k)
        # Get a vector of euclidean distances
        distance <- sapply(seq_len(nrow(cells)), function(i) {
            dist(rbind(cells[,test.markers][i,], imputed.cells[i,]))
        })
    })

    # Get final distances into tibble format
    names(final.distances) <- k.titration
    final.distances <- as.tibble(final.distances)
    final.distances <- apply(final.distances, 2, mean)

    return(final.distances)
}
