#' @import rflann
NULL

#' @title Compute knn using the fnn algorithm
#'
#' @description This function is a wrapper around FNN package
#' functionality to speed up the KNN process. It uses KD trees as default,
#' along with k set to 100. Selection of K will vary based on the dataset.
#' See k.selection.R.
#' @param cell.df the cell data frame used as input
#' @param input.markers markers to be used as input for knn
#' @param k the number of nearest neighbors to identify
#' @return nn: list of 2, nn.index: index of knn (columns) for each cell (rows)
#' nn.dist: euclidean distance of each k-nearest neighbor
#' @examples
#' fnn(combined[1:1000,], input.markers)
#' @export
fnn <- function(cell.df, input.markers, k = 100) {
    print("finding k-nearest neighbors")

    # Edge case (rflann with kd-tree doesn't have it)
    if(k >= nrow(cell.df)) {
        stop("k must be less than the total number of data points")
    }

    input <- cell.df[,input.markers]

    # Using the rflann package
    nn <- Neighbour(query = input, ref = input, k = k + 1)
    nn.index <- nn[[1]][,2:ncol(nn[[1]])]
    nn.dist <- nn[[2]][,2:ncol(nn[[2]])]

    print("k-nearest neighbors complete")
    return(list(nn.index = nn.index, nn.dist = nn.dist))
}


#' @title Performs a series of statistical tests on the batch of cells
#' of interest.
#'
#' @description This function performs the statistics across the nearest
#' neighborhoods, and is one of the main workhorses within the scone.values
#' function
#'
#' @param basal tibble of cells corresponding to the unstimulated condition
#' @param stim a tibble of cells corresponding to the stimulated condition
#' @param fold a string that specifies the use of "median" or "mean" when
#' calculating fold change
#' @param stat.test a string that specifies Mann-Whitney U test (mwu) or T test (t)
#' for q value calculation
#' @param stim.name a string corresponding to the name of the stim being tested
#' compared to basal
#' @return result: a named vector corresponding to the results of the
#' "fold change" and mann-whitney u test
run.statistics <- function(basal,
                           stim,
                           fold = "median",
                           stat.test = "mwu",
                           stim.name) {
    # Edge case of a knn consisting of only one of the two conditions
    # More common in messier datasets
    if(nrow(basal) == 0 | nrow(stim) == 0) {
        return(rep(NA, times = 2*ncol(basal) + 1))
    }

    # Fold change between unstim and stim (output is named vector)
    if(fold == "median") {
        fold <- apply(stim, 2, median) - apply(basal, 2, median)
    } else if (fold == "mean") {
        fold <- apply(stim, 2, mean) - apply(basal, 2, mean)
    } else {
        stop("please select median or mean to be used as input for raw
             change calculation")
    }

    # Mann-Whitney U test or T test
    if(stat.test == "mwu") {
        qvalue <- sapply(1:ncol(basal), function(j) {
            p <- wilcox.test(basal[[j]], stim[[j]])$p.value
            return(p)
        })
    } else if (stat.test == 't') {
        qvalue <- sapply(1:ncol(basal), function(j) {
            p <- t.test(basal[[j]], stim[[j]])$p.value
            return(p)
        })
    } else {
        stop("please select either Mann-Whitney U test (mwu) or T test (t) for input")
        return()
    }

    # Naming the vectors
    names(qvalue) <- names(fold) # qvalue is not yet a named vector, so its named here
    names(qvalue) <- paste(names(qvalue), stim.name, "qvalue", sep = ".") # specifying
    names(fold) <- paste(names(fold), stim.name, "change", sep = ".") # specifying

    # Get the unstim and stim thresholds done
    fraction.cond2 <- nrow(stim)/sum(nrow(basal), nrow(stim))
    names(fraction.cond2) <- paste(stim.name, "fraction.cond.2", sep = ".")
    result <- c(qvalue, fold, fraction.cond2)
    return(result)
}


#' @title Runs a t test on the medians or means of multiple donors for the
#' same condition
#'
#' @description This function is for the instance that multiple donors are
#' being compared against each other within the k-nearest neighborhood of
#' interest. The mean value of the markers of interest is calculated across
#' the donors, such that each data point for the subsequent t-test represents
#' a marker for a danor.
#' @param basal tibble that contains unstim for a knn including donor identity
#' @param stim tibble that contains stim for a knn including donor identity
#' @param stim.name string of the name of the current stim being tested
#' @param donors vector of strings corresponding to the designated names
#' of the donors
#' @return result: a named vector of p values (soon to be q values) from the
#' t test done on each marker
multiple.donor.statistics <- function(basal, stim, stim.name, donors) {
    # get the means of all the donors

    basal.stats <- tibble()
    stim.stats <- tibble()
    for(i in donors) {
        basal.curr <- basal[basal$donor == i,] %>%
            .[!(colnames(.) %in% "donor")]
        stim.curr <- stim[stim$donor == i,] %>%
            .[!(colnames(.) %in% "donor")]

        basal.mean <- apply(basal.curr, 2, mean)
        stim.mean <- apply(stim.curr, 2, mean)

        basal.stats <- rbind(basal.stats, basal.mean)
        stim.stats <- rbind(stim.stats, stim.mean)
    }
    # assumes "donor" always placed at end!
    colnames(basal.stats) <- colnames(basal)[-ncol(basal)]
    colnames(stim.stats) <- colnames(stim)[-ncol(basal)]

    # T testing (only if there's no NA in the dataset)
    if(nrow(na.omit(basal.stats)) == length(donors) & # Watch the newline here
       nrow(na.omit(stim.stats)) == length(donors)) {
        result <- sapply(1:ncol(basal.stats), function(i) {
            t.test(basal.stats[[i]], stim.stats[[i]])$p.value
        })
    } else {
        result <- rep(NA, times = ncol(basal.stats))
    }


    names(result) <- paste(colnames(basal.stats),
                           stim.name,
                           "replicate.qvalue",
                           sep = ".")
    return(result)
}

#' @title Corrects all p values for multiple hypotheses, sets threshold for
#' which change values should be reported
#'
#' @description Given the number of comparisons we make across k-nearest
#' neighborhoods, which is far more than that of disjoint subsetting, this
#' step is important given that there is an increased likelihood that some
#' statistically significant differences will occur by chance.
#' @param cells tibble of change values, p values, and fraction condition 2
#' @param threshold a q value below which the change values will be reported
#' for that cell for that param. If no change is desired, this is set to 1.
#' @return inputted p values, adjusted and therefore described as "q values"
q.correction.thresholding <- function(cells, threshold) {
    # Break apart the result
    fold <- cells[,grep("change$", colnames(cells))]
    qvalues <- cells[,grep("qvalue$", colnames(cells))]
    ratio <- cells[,grep("cond2$", colnames(cells))]

    # rest <- cells[,!(colnames(cells) %in% colnames(qvalues))]

    # P value correction
    qvalues <- apply(qvalues, 2, function(x) p.adjust(x, method = "BH")) %>%
        as.tibble

    # Thresholding the raw change
    if(threshold < 1) {
        names <- colnames(fold)
        fold <- lapply(1:ncol(fold), function(i) {
            curr <- fold[[i]]
            curr <- ifelse(qvalues[[i]] < threshold, curr, 0)
        }) %>% do.call(cbind, .) %>%
            as.tibble()
        colnames(fold) <- names
    }

    #Bring it all together
    result <- bind_cols(qvalues, fold, ratio)
    return(result)
}

#' @title Get the KNN density estimaion
#' @description Obtain a density estimation derived from the original manifold,
#' avoiding the lossiness of lower dimensional embeddings
#' @param nn.matrix A list of 2, where the first is a matrix of nn indices,
#' and the second is a matrix of nn distances
#' @return a vector where each element is the KNN-DE for that given cell,
#' ordered by row number, in the original input matrix of cells x features
#' @examples
#' ex.knn <- fnn(combined[1:1000,], input.markers)
#' get.knn.de(ex.knn)
#' @export
get.knn.de <- function(nn.matrix) {
    # Distances
    nn.dist <- nn.matrix[[2]]

    # KNN-DE
    mean.dist <- apply(nn.dist, 1, mean)
    density <- 1/mean.dist
    return(density)
}

#' @title Make list of cells by features for each KNN member
#' @description Takes the KNN function output and the cell data, and
#' makes list where each element is a matrix of cells in the KNN and features.
#' @param cell.data tibble of cells by features
#' @param nn.matrix list of 2. First element is cells x 100 nearest neighbor
#' indices. Second element is cells x 100 nearest neighbor distances
#' @return a list where each element is the cell number from the
#' original cell.data tibble and a matrix of cells x feautures for its KNN
#' @examples
#' ex.knn <- fnn(combined[1:1000,], input.markers)
#' make.knn.list(combined[1:1000,], ex.knn)
#' @export
make.knn.list <- function(cell.data, nn.matrix) {
    # Unpack the KNN output
    nn.index <- nn.matrix[[1]]

    # The list
    knn.list <- lapply(1:nrow(nn.index), function(i) {
        cell.data[nn.index[i,],]
    })

    return(knn.list)
}


#' @title Master function for per-knn statistics functionality, integrating the
#' other non-exported functions within this script.
#'
#' @description This function is run following the KNN computation
#' and respective cell grouping. The function also contains a progress ticker
#' that allows one to determine how much time left in this function.
#' @param nn.matrix a matrix of cell index by nearest neighbor index, with
#' values being cell index of the nth nearest neighbor
#' @param cell.data tibble of cells by features
#' @param scone.markers vector of all markers to be interrogated via
#' statistical testing
#' @param unstim an object (used so far: string, number)
#' specifying the "basal" condition
#' @param threshold a number indicating the p value the raw change should
#' be thresholded by.
#' @param fold a string that specifies the use of "median" or "mean" when
#' calculating fold change
#' @param stat.test string denoting Mann Whitney U test ("mwu") or T test ("t)
#' @param multiple.donor.compare a boolean that indicates whether t test
#' across multiple donors should be done
#' @return result: tibble of raw changes and p values for each feature of
#' interest, and fraction of cells with condition 2
#' @examples
#' ex.cells <- combined[sample(1:nrow(combined), 1000),]
#' ex.nn <- fnn(ex.cells, input.markers)
#' scone.values(ex.nn, ex.cells, funct.markers, "basal")
#' @export
scone.values <- function(nn.matrix,
                         cell.data,
                         scone.markers,
                         unstim,
                         threshold = 0.05,
                         fold = "median",
                         stat.test = "mwu",
                         multiple.donor.compare = FALSE) {

    print("running per-knn statistics")

    # Get the donor names if you need it
    if(multiple.donor.compare == TRUE) {
        donors <- unique(cell.data$donor)
    }

    # Get all stim.names
    conditions <- unique(cell.data$condition)
    stims <- conditions[!(conditions %in% unstim)]

    # Unpack the nn matrix
    nn.dist <- nn.matrix[[2]]
    nn.matrix <- nn.matrix[[1]] # overwrite to be used in rest of this

    # Get the KNN density estimation (1 / avg distance to knn)
    avg.dist <- apply(nn.dist, 1, mean)
    knn.density <- 1/avg.dist

    # Variables to be used for progress tracker within the loop
    percent <- nrow(cell.data) %/% 10

    final <- lapply(stims, function(s) {
        # Process each column in the nn matrix
        # This makes a list of vecors corresponding qvalue and fold change
        count <- 0
        result <- lapply(1:nrow(nn.matrix), function(i) {
            # A tracker
            if(i %% percent == 0) {
                count <<- count + 10
                print(paste(count, "percent complete", sep = " "))
            }

            # Index the nn matrix
            curr <- cell.data[nn.matrix[i,],]
            basal <- curr[curr$condition == unstim,] %>% .[,scone.markers]
            stim <- curr[curr$condition == s,] %>% .[,scone.markers] # Change this to specify stim name

            # Fold change cmoparison and Mann-Whitney U test, along with "fraction condition 2"
            output <- run.statistics(basal, stim, fold, stat.test, s)

            # Note that this will overwrite the initial output (for now)
            if(multiple.donor.compare == TRUE) {
                nn.donors <- unique(curr$donor)

                # We want multiple donor testing only if all donors are in each knn
                if(length(unique(nn.donors)) < length(donors)) {
                    donor.output <- rep(NA, times = length(scone.markers)) # Check length
                } else {
                    basal.d <- curr[curr$condition == unstim,] %>% .[,c(scone.markers, "donor")]
                    stim.d <- curr[curr$condition == s,] %>% .[,c(scone.markers, "donor")]
                    donor.output <- multiple.donor.statistics(basal.d, stim.d, s, donors)
                }
                #names(donor.output) <- paste(scone, s, "donor.t.test.qvalue", sep = ".") # Name the donor t test vector
                output <- c(output, donor.output)
            }

            return(output)
        })

        # Melt the list together into a tibble and return it
        result <- do.call(rbind, result) %>%
            as.tibble()

        # print(result[,grep("Stat3", colnames(result))])
        return(result)
    })

    # Returning an error message for the Zunder dataset, but not the Wanderlust dataset
    final <- do.call(cbind, final) %>%
        as.tibble()

    # Do the p value correction, followed by thresholding if you're going to
    frac <- final[, grep("fraction", colnames(final))]
    final <- q.correction.thresholding(final, threshold)
    final <- bind_cols(final, frac)

    # Add the density estimation
    final$density <- knn.density

    print("per-knn statistics complete")

    # Change the only "character" column to a numeric column
    return(final)
}
