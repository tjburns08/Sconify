#' @import Rtsne ggplot2
#' @title Add tSNE to your results.
#'
#' @description This function gives the user the option to add t-SNE to the
#' final output, using the same input features used in KNN (eg. surface markers)
#' as input for t-SNE.
#' @param dat: matrix of cells by features, that contain all features needed
#' for tSNE analysis
#' @param input: the features to be used as input for tSNE (usually the same
#' for knn generation)
#' @return result: dat, with tSNE1 and tSNE2 attached
add.tsne <- function(dat, input) {
    result <- Rtsne(X = dat[,input],
                    dims = 2,
                    pca = FALSE,
                    verbose = TRUE)$Y %>%
        as.tibble
    names(result) <- c("bh-SNE1", "bh-SNE2")
    result <- bind_cols(dat, result)
    return(result)
}


#' @title Log transform the q values
#'
#' @description Takes all p values from the data and does a log10 transform
#' for easier visualization.
#' @param dat: tibble containing cells x features, with orignal expression,
#' p values, and raw change
#' @param negative: boolean value to determine whether to multiple transformed
#' p values by -1
#' @return result: tibble of cells x features with all p values log10
#' transformed
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

#' @title Transform strings to numbers.
#'
#' @description Takes a vector of strings and outputs simple numbers. This
#' takes care of the case where conditions are listed as strings (basal, IL7),
#' in which case they are converted to numbers (1, 2)
#' @param strings vector of strings
#' @return strings: same vector with each unique element converted to a number
string.to.numbers <- function(strings) {
    elements <- unique(strings)
    for(i in 1:length(elements)) {
        strings <- ifelse(strings == elements[i], i, strings)
    }
    return(as.numeric(strings))
}


#' @title Post-processing for scone analysks.
#'
#' @description Performs final processing and transformations on the scone data
#' @export
#' @param scone.output: tibble of the output of the given scone analysis
#' @param cells: the tibble used as input for the scone.values function
#' @param input: the input markers used for the knn calculation (to be used
#' for tsne here)
#' @param tsne: boolean value to indicate whether tSNE is to be done
#' @param log.transform.qvalue: boolean to indicate whether log transformation
#' of all q values is to be done
#' @return result: the concatenated original input data with the scone derived
#' data, with the option of the q values being inverse log10 transformed, and
#' two additional tSNE columns being added to the data (from the Rtsne package)
#' @examples
#' final <- post.processing(scone.output = scone.output,
#'                          cell.data = combined,
#'                          input = input.markers)
#' @export
post.processing <- function(scone.output,
                            cell.data,
                            input,
                            tsne = TRUE,
                            log.transform.qvalue = TRUE) {
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

#' @title Takes as input a labeled data matrix,
#' and outputs an fcs file containing said inroamtion.
#'
#' @description Note that Cytobank's DROP feature in theory does what
#' this script does for csv files, but its availability is currently limited
#' and not all fcs processing software will have similar functionality.
#' Note also the untransform parameter. This is for the instrances where
#' an asinh transformation with a scale argument of 5 takes place within R
#' (for example when one is clustering). Cytobank is not used to reading these
#' scales. Thus, untransforming them prior to Cytobank use allows the biaxials
#' to be effectively visualized.
#'
#' @param dat a data matrix intended to be converted to fcs
#' @param outfile a string containing your original data matrix
#' @param untransform takes data out of asinh(x/5) transformation
#' @return fcs file containing your original data matrix
#' @examples
#' basal <- system.file('extdata',
#'     'Bendall et al Cell Sample C_basal.fcs',
#'     package = "Sconify")
#' basal <- fcs.to.tibble(basal, transform = "asinh")
#' data.to.fcs(basal, "basal.output.FCS")
#' @export
data.to.fcs <- function(dta, outfile, untransform) {
    # Convert the data of interest into a matrix
    if(!is.matrix(dta)) {
        dta <- as.matrix(dta)
    }

    # Take it out of asinh transform
    if(untransform == TRUE) {
        dta <- sinh(dta)*5
    }

    # Convert into a flow frame
    dta <- flowFrame(dta)

    # now you can save it back to the filesystem
    write.FCS(dta, outfile)
}

#' @title make.hist
#'
#' @description Makes a histogram of the data that is inputted
#'
#' @param dat tibble consisting both of original markers and the appended values from scone
#' @param k the binwidth, set to 1/k
#' @param column.label the label in the tibble's columns the function will search for
#' @param x.label the label that the x axis will be labeled as
#'
#' @return a histogram of said vector in ggplot2 form
#' @export
make.hist <- function(dat,
                      k,
                      column.label,
                      x.label) {
    ggplot(data = dat, aes(x = dat[,grep(column.label, colnames(dat))])) +
        geom_histogram(aes(y = ..count..), binwidth = 1/k) +
        xlim(c(0, 1)) +
        theme(text = element_text(size = 20)) +
        xlab(x.label)
}
