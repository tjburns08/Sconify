#' @import magrittr tibble flowCore
NULL

#' @title Takes in an example file as input and returns all markers
#'
#' @description This is a quick way to get a list of preferred marker names.
#' This outputs a csv file containing all markers in the dataset in the name
#' format that will be recognized by downstream functions. You manually
#' alter this list to remove and/or categorize the said markers. The file
#' can then be read in (stringsAsFactors = FALSE) to give you the marker
#' list of interest. In particular, name the top of the column as "markers" if
#' you're just altering the list. If you're doing to divide it into static
#' and functional markers, produce two columns, naming them respectively.
#'
#' @param file the fcs file of interest
#' @return the list of markers of interest, as a csv
#' @examples
#' file <- system.file("extdata", "Bendall_et_al_Cell_Sample_C_basal.fcs",
#'     package = "Sconify")
#' get.marker.names("file")
#' @export
get.marker.names <- function(file) {
    cells <- fcs.to.tibble(file)
    result <- colnames(cells)
    write.csv(result, "markers.csv", row.names = FALSE)
    return(result)
}

#' @title Takes a file as input and returns a data frame of cells by features
#' @description This function is a quick way to take the exprs content of
#' a fcs file, do an asinh transform, and create a tibble data structure
#' that can be further manipulated. Our default transform is asinh, but
#' you just have to change the transform to anything else, and you'll
#' get the raw data. This function is used in the main function
#' process.multiple.files
#'
#' @param file the fcs file containing cell infomration
#' @param transform if set to asinh, then asinh transforms with scale arg 5
#' @return tibble of info contained within the fcs file
#' @examples
#' file <- system.file("extdata", "Bendall_et_al_Cell_Sample_C_basal.fcs",
#'     package = "Sconify")
#' fcs.to.tibble(file)
#' @export
fcs.to.tibble <- function(file, transform = "asinh") {
    # Read in the files and set the columns as human-named parameters
    cells <- flowCore::read.FCS(file = file)
    params <- as.vector(pData(parameters(cells))$desc)
    colnames(cells) <- params

    # Turn into matrix and asinh tranform with cofactor of 5
    if(transform == "asinh") {
        cells <- exprs(cells) %>%
            `/`(5) %>%
            asinh() %>%
            as_tibble()
    } else {
        cells <- exprs(cells) %>%
            as_tibble()
    }
    return(cells)
}

#' @title Performs quantile normalization on the data frame (patient) of interest
#'
#' @description Credit goes to: http://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#'
#' @param df a data frame with rows as cells and columns as features
#' @return a data frame where the columns have been quantile normalized
quantile_normalization <- function(df){
    df_rank <- apply(df,2,rank,ties.method="first")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)

    index_to_mean <- function(my_index, my_mean){
        return(my_mean[my_index])
    }

    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(as.data.frame(df_final))
}

#' @title Takes a list of tibbles as input, and performs per-column
#' quantile normalization, then outputs the quantile normalized list
#' @description This function performs per-marker quantile normalization
#' on multiple data tibbles. The normalization occurrs marker by marker.
#' The user assumes that the markers are distributed equally across tibbles,
#' as quantile normalization forces these marker distributions to be the same
#' per file
#'
#' @param dat.list a list of tibbles
#' @return the per-column quantile normalized list
quant.normalize.elements <- function(dat.list) {
    # Store the marker names for the re-naming when the list is reverted
    marker.names <- colnames(dat.list[[1]])

    # Transpose the list such that it's now marker-by-marker,
    # and quantile normalize
    new.list <- lapply(1:ncol(dat.list[[1]]), function(i){
        curr <- lapply(dat.list, function(j) {
            slice <- j[,i]
        })
        curr <- do.call(what = cbind, args = curr) %>% quantile_normalization()
    })

    # Transpose again, leading to the original list of tibbles
    revert <- lapply(1:ncol(new.list[[1]]), function(i){
        curr <- lapply(new.list, function(j) {
            slice <- j[,i]
        })
        curr <- do.call(what = cbind, args = curr) %>% as.tibble()
    })

    # Loop through each element of the list and re-name the columns accordingly
    for(i in 1:length(revert)) {
        colnames(revert[[i]]) <- marker.names
    }

    return(revert)
}


#' @title Converts multiple files into a concatenated data frame
#'
#' @description This is tailored to a very specific file format for unstim/stim
#' Files need the following name convention: "xxxx_stim.fcs"
#' Files where you want to name the patients need the following convention:
#' "xxxx__patientID_stim.fcs"
#' @param files a vector of file names (name = "anything_condition.fcs")
#' @param transform set to asinh if you want to do an asinh transform of
#' all markers in the dataset
#' @param numcells desiered number of cells in the matrix, set at 10k
#' @param norm boolean that quantile normalizes the data if true
#' @param scale boolean that converts all values to z scores if true
#' @param input the static markers that will be used downstream in knn
#' computation. These are included here to include the option of per-marker
#' quantile normalization, in the event norm is set to TRUE
#' @param name.multiple.donors boolean indicating whether multiple donors
#' will be distinguished (as a separate "patient" column)
#' @return result: a combined file set
#' @examples
#' file1 <- system.file("extdata", "Bendall_et_al_Cell_Sample_C_basal.fcs",
#'     package = "Sconify")
#' file2 <- system.file("extdata", "Bendall_et_al_Cell_Sample_C_IL7.fcs",
#'     package = "Sconify")
#' process.multiple.files(c(file1, file2), input.markers)
#' @export
process.multiple.files <- function(files,
                                   transform = "asinh",
                                   numcells = 10000,
                                   norm = FALSE,
                                   scale = FALSE,
                                   input,
                                   name.multiple.donors = FALSE) {
    n <- numcells %/% length(files) # integer division

    # Turn input into a list
    dat <- lapply(files, function(i){
        curr <- fcs.to.tibble(i, transform = transform) %>%
            .[,na.omit(colnames(.))] %>%
            as.tibble() # in case columns are named "NA"

        # Subsample if the number of cells is greater than specified n
        if(nrow(curr) > n) {
            curr <- curr[sample(nrow(curr), n),]
        }

        # Break curr into the values you're going to normalize for
        # knn generation and those you're not
        curr.input <- curr[,input]
        curr.rest <- curr[,(!(colnames(curr) %in% input))]


        if(scale == TRUE) {
            curr.input <- apply(curr.input, 2, scale) %>% as.tibble()
            # takes care of scaling where all values are same
            curr.input <- replace(curr.input, is.na(curr.input), 1)
        }

        # Overwrite curr, bringing the (potentially normalized/scaled)
        # curr.input and curr.rest together
        curr <- bind_cols(curr.input, curr.rest)

        # File standard: "anything_condition.fcs"
        curr$condition <- sub(".*_", "", i) %>% sub("\\..*", "", .)

        if(name.multiple.donors == TRUE) {
            curr$donor <- sub(".*__", "", i) %>% sub("\\_.*", "", .)
        }

        return(curr)
    })

    if(norm == TRUE) {
        if(length(files) < 2) {
            stop("Quantile normalization can only happen with 2 or more files")
        }
        # The list with only the input markers,
        # and perform quantile normalization
        dat.input <- lapply(dat, function(i) {
            i[,input]
        }) %>% quant.normalize.elements(.)

        # The list with only the non-input markers
        dat.rest <- lapply(dat, function(i) {
            i[,(!(colnames(i) %in% input))]
        })

        # Merge the two lists together again
        dat <- lapply(1:length(dat), function(i) {
            bind_cols(dat.input[[i]], dat.rest[[i]])
        })
    }

    # Bind the list of tibbles by the row
    if(length(dat) > 1) {
        result <- bind_rows(dat)
    } else {
        result <- dat[[1]]
    }

    return(result)
}


#' @title Runs "process.multiple.files" on a single file, splits it randomly,
#' and presends half of it is "unstim" and half of it is "stim"
#'
#' @description This is meant to serve as a control for the basic "unstim"
#' and "stim" pipeline that is generally used within this package.
#' If phosphoproteins are being compared across conditions, for example,
#' then there should be no difference in the case that you split the same file
#' and compare the two halves.
#'
#' @param file the file we're going to split
#' @param transform if set to asinh, performs asinh transformation on all
#' markers of the dataset
#' @param numcells the number of total cells to be subsampled to, set at 10k
#' for default
#' @param norm boolean of whether data should be quantile normalized
#' @param scale boolean of whether data should be z-scored
#' @param input.markers vector of strings indicating the markers
#' to be used as input
#' @return tibble containing original markers and all values
#' calculated by SCONE
#' @examples
#' file <- system.file("extdata", "Bendall_et_al_Cell_Sample_C_basal.fcs",
#'     package = "Sconify")
#' splitFile(file, input.markers)
#' @export
splitFile <- function(file,
                      transform = "asinh",
                      numcells = 10000,
                      norm = FALSE,
                      scale = FALSE,
                      input.markers) {
    # Create a subsample
    total.unstim <- process.multiple.files(files = file,
                                           numcells = numcells,
                                           transform = transform,
                                           norm = norm,
                                           scale = scale,
                                           input = input.markers,
                                           name.multiple.donors = FALSE)

    cell.rows <- 1:nrow(total.unstim)

    # If the number of cells is odd, get rid of the bottom cell in the matrix
    if(nrow(total.unstim) %% 2 != 0) {
        cell.rows <- cell.rows[-length(cell.rows)]
    }

    # Random sampling
    index <- sample(cell.rows, size = nrow(total.unstim)/2, replace = FALSE)
    index.2 <- cell.rows[!(cell.rows %in% index)]
    s1 <- total.unstim[index,]
    s1$condition <- "Split.1"
    s2 <- total.unstim[index.2,]
    s2$condition <- "Split.2"
    s.merge <- bind_rows(s1, s2)
    return(s.merge)
}

#' @title Meaning of life
#'
#' @description Just a random musing
#'
#' @examples
#' meaning.of.life()
#' @export
meaning.of.life <- function() {
    exist
}

