############################### SELECTION OF IDEAL K ###############################

# Imputes values for all markers (used as input) for each cell
# Args:
#   cells: the input matrix of cells
#   input.markers: the markers the user wants to impute
#   nn: the matrix of k-nearest neighbors (derived perhaps NOT from the "input markers" above)
# Returns:
#   result: a data frame of imputed cells for the "input markers" of interest
impute <- function(cells, input.markers, nn) {
    result <- sapply(1:nrow(cells), function(i) {
        curr.nn <- cells[nn[i,],][,input.markers] %>% apply(., 2, median)
    }) %>%
        t() %>%
        as.data.frame()
    return(result)
}

# Tests the euclidean distance error for imputation using knn and markers of interest
# Args:
#   k.titration: a vector integer values of k to be tested
#   cells: a matrix of cells by features used as original input
#   input.markers: markers to be used for the knn calculation
#   test.markers: the markers to be tested for imputation (either surface or scone)
# Returns:
#   final.distances: the median imputation error for each value k tested
impute.testing <- function(k.titration, cells, input.markers, test.markers) {
    final.distances <- lapply(k.titration, function(k) {
        nn <- fnn(cell.df = cells, input.markers = input.markers, k = k)
        imputed.cells <- impute(cells = cells, input.markers = test.markers, nn = nn)

        print(k)
        # Get a vector of euclidean distances
        distance <- sapply(1:nrow(cells), function(i) {
            dist(rbind(cells[,test.markers][i,], imputed.cells[i,]))
        })
    })

    # Get final distances into tibble format
    names(final.distances) <- k.titration
    final.distances <- as.tibble(final.distances) # BUG
    final.distances <- apply(final.distances, 2, median)

    return(final.distances)
}
