############################### POST-PROCESSING AND TSNE ###############################

# Adds tSNE1 and tSNE2 to the columns of a dataset
# Args:
#   dat: matrix of cells by features, that contain all features needed for tSNE analysis
#   input: the features to be used as input for tSNE (usually the same for knn generation)
# Returns:
#   result: dat, with tSNE1 and tSNE2 attached
add.tsne <- function(dat, input) {
    result <- Rtsne(X = dat[,input], dims = 2, pca = FALSE, verbose = TRUE)$Y %>% as.tibble
    names(result) <- c("bh-SNE1", "bh-SNE2")
    result <- bind_cols(dat, result)
    return(result)
}


# Takes all p values from the data and does a log10 transform for visualization
# Args:
#   dat: tibble containing cells x features, with orignal expression, p values, and raw change
#   negative: boolean value to determine whether to multiple transformed p values by -1
# Returns:
#   result: tibble of cells x features with all p values log10 transformed
log.transform.q <- function(dat, negative) {

    # Split the input
    qvalue <- dat[,grep("qvalue$", colnames(dat))]
    rest <- dat[, !(colnames(dat) %in% colnames(qvalue))]
    qvalue <- apply(qvalue, 2, log10) %>% as.tibble()

    # Inverse sometimes better for visualization
    if(negative == TRUE) {
        qvalue <- apply(qvalue, 2, function(x) x*-1) %>% as.tibble()
    }

    # Note that the order will be slightly different that the original "final"
    result <- bind_cols(rest, qvalue)
    return(result)
}

# Takes a vector of strings and outputs simple numbers
# Args:
#   strings: vector of strings
# Returns:
#   strings: same vector with each unique element converted to a number
string.to.numbers <- function(strings) {
    elements <- unique(strings)
    for(i in 1:length(elements)) {
        strings <- ifelse(strings == elements[i], i, strings)
    }
    return(as.numeric(strings))
}


# Performs final processing and transformations on the scone data
# Args:
#   scone.output: tibble of the output of the given scone analysis
#   cells: the tibble used as input for the scone.values function
#   input: the input markers used for the knn calculation (to be used for tsne here)
#   tsne: boolean value to indicate whether tSNE is to be done
#   log.transform.qvalue: boolean to indicate whether log transformation of all q values is to be done
# Returns:
#   result: the concatenated original input data with the scone derived data, with the option of
#     the q values being inverse log10 transformed, and two additional tSNE columns being added
#     to the data (from the Rtsne package)
post.processing <- function(scone.output, cell.data, input, tsne = TRUE, log.transform.qvalue = TRUE) {
    # Generic pre-processing
    result <- bind_cols(cell.data, scone.output) %>% na.omit()
    result$condition <- string.to.numbers(result$condition)

    # Adding two tSNE columns
    if(tsne == TRUE) {
        result <- add.tsne(dat = result, input = input)
    }

    # Doing an inverse log transformation of the q value
    if(log.transform.qvalue == TRUE) {
        result <- log.transform.q(dat = result, negative = TRUE)
    }

    return(result)
}

