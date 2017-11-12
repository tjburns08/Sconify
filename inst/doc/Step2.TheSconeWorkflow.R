## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE, warning = FALSE)

## ------------------------------------------------------------------------
# The above example is stored in the package as "markers," loaded here
library(Sconify)
library(dplyr)

# How to convert your excel sheet into vector of static and functional markers
data(markers)

# addition in the end gets rid of empty strings
input.markers <- markers$surface %>% .[. != ""] 
funct.markers <- markers$functional %>% .[. != ""]

# Selection of the k. Default is number of cells / 100. 
k <- nrow(combined) %/% 100

# The built-in scone functions
nn <- fnn(cell.df = combined, input.markers = input.markers, k = k)

# Cell identity is in rows, k-nearest neighbors are columns
# List of 2 includes the cell identity of each nn, and the euclidean distance between
#   itself and the cell of interest
str(nn)
nn[1:20, 1:10]

## ------------------------------------------------------------------------
scone.output <- scone.values(nn.matrix = nn, 
                      cell.data = combined, 
                      scone.markers = funct.markers, 
                      unstim = "basal")
scone.output

