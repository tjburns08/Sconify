# Sconify
### Continuous visualization of differences between biological conditions in single-cell data

In high-dimensional single cell data, comparing changes in functional markers between conditions is typically done across manual or algorithm-derived partitions based on population-defining markers. This package performs these comparisons across overlapping k-nearest neighbor (KNN) groupings. Each cell thus represents a proxy of its local neighborhood, and is colored by a comparison of interest for visualization in low-dimensional embeddings (eg. t-SNE). This package includes an objective optimization of k based on minimizing functional marker KNN imputation error. Proof-of-concept work has visualized the exact location of an IL-7 responsive subset in a B cell developmental trajectory on a t-SNE map independent of clustering. Cell frequency analysis revealed that KNN is sensitive to detecting artifacts due to marker shift, and therefore can also be valuable in one’s quality control pipeline. Overall, we found that KNN groupings can efficiently extract a large amount of information from mass cytometry data, making it useful especially for the initial stages of data analysis. 

## Installation

The pacakge was recently accepted into BioConductor, so this section will contain new instructions shortly (though the content of the package will remain the same). For now, please install the package directly through GitHub as follows:

```
library(devtools)
devtools::install_github("tjburns08/Sconify", build_vignettes = TRUE)
```

## Further instruction

This package contains multiple vignettes, explaining all aspects of the Sconify package, with plenty of examples and pictures. Once the package has been installed, please type the following:

```
library(Sconify)
browseVignettes("Sconify")
```

