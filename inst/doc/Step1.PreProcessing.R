## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE, warning = FALSE)

## ---- out.width = "200px", eval = TRUE-----------------------------------
knitr::include_graphics("original.markers.csv.example.png")


## ---- out.width = "200px"------------------------------------------------
knitr::include_graphics("modified.markers.csv.example.png")


## ------------------------------------------------------------------------

library(Sconify)

# FCS file provided in the package
basal <- system.file('extdata',
    'Bendall et al Cell Sample C_basal.fcs',
    package = "Sconify")

# Example of data with no transformation
basal.raw <- fcs.to.tibble(basal, transform = "none")
basal.raw

# Asinh transformation with a scale argument of 5
basal.asinh <- fcs.to.tibble(basal, transform = "asinh")
basal.asinh


## ------------------------------------------------------------------------
# The FCS files (THEY NEED TO BE IN THE "....._condidtion.fcs" format")
basal <- system.file('extdata',
    'Bendall et al Cell Sample C_basal.fcs',
     package = "Sconify")
il7 <- system.file('extdata',
    'Bendall et al Cell Sample C_IL7.fcs',
    package = "Sconify")

# The markers (after they were modified by hand as instructed above)
markers <- system.file('extdata',
    'markers.csv',
    package = "Sconify")
markers <- read.csv(markers, stringsAsFactors = FALSE)
surface <- markers$surface

# Combining these. Note that default is sub-sampling to 10,000 cells, not normalizing, and not scaling
combined <- process.multiple.files(files = c(basal, il7), input = surface)
combined
unique(combined$condition)

# Limit your matrix to surface markers, if just using those downstream
combined.input <- combined[,surface]
combined.input

# We can do this on a single file as well. Notice I allow for scaling and I change the number of cells
basal.data <- process.multiple.files(files = basal, numcells = 6000, scale = TRUE, input = surface)
basal.data
unique(basal.data$condition)


## ------------------------------------------------------------------------
# The FCS files
basal <- system.file('extdata',
    'Bendall et al Cell Sample C_basal.fcs',
    package = "Sconify")
markers <- system.file('extdata',
    'markers.csv',
    package = "Sconify")

# The markers
markers <- read.csv(markers, stringsAsFactors = FALSE)
surface <- markers$surface

# The function
split <- splitFile(file = basal, input.markers = surface)
split
unique(split$condition)


