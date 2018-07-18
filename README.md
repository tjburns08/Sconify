# Sconify
### Continuous visualization of differences between biological conditions in single-cell data

In high-dimensional single cell data, comparing changes in functional markers between conditions is typically done across manual or algorithm-derived partitions based on population-defining markers. This package performs these comparisons across overlapping k-nearest neighbor (KNN) groupings. Each cell thus represents a proxy of its local neighborhood, and is colored by a comparison of interest for visualization in low-dimensional embeddings (eg. t-SNE). This package includes an objective optimization of k based on minimizing functional marker KNN imputation error. Proof-of-concept work has visualized the exact location of an IL-7 responsive subset in a B cell developmental trajectory on a t-SNE map independent of clustering. Cell frequency analysis revealed that KNN is sensitive to detecting artifacts due to marker shift, and therefore can also be valuable in one’s quality control pipeline. Overall, we found that KNN groupings can efficiently extract a large amount of information from mass cytometry data, making it useful especially for the initial stages of data analysis. 

## Installation

The Sconify package is available on BioConductor. You can install it as 
follows:

```
# Install the BioConductor manager from CRAN
install.packages("BiocManager")

# Install the Sconify package from BioConductor
BiocManager::install("Sconify")
```

The development version of Sconify can be found on Github. You can install
it as follows:

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

## Visual explanations
![Alt text](vignettes/sconify_visual_explanation.png?raw=true "Title")

This is the general schematic of Sconify. Concatenated data in high dimensional space is grouped with each cell's k-nearest neighbors. Statistics are performed within each neighborhood. This allows for dimension reduction maps, like t-SNE to have the "fold change functionality" that has been desired since t-SNE became a staple of CyTOF analysis pipelines. 

![Alt text](vignettes/sconify_proof_of_concept.png?raw=true "Title")

This is an example of Sconify in use. Notice on the left there is an "untreated" and "IL7" t-SNE map, followed by a Sconify-enabled composite t-SNE map that shows the fold-change values. The Sconify package can visualize fold change, p-values for t and Mann-Whitney U tests, and pvalue-thresholded fold change. Note that t-SNE can be run directly within Sconify, and the output can be read in, as a csv, to Cytobank or CYT for further visualization (aside from the ggplot-based visualization that Sconify provides. 

## For those of you with little programming experience

There is a lite version of Sconify under construction at www.sconify.org. It can currently handle small fcs files, and it will be updated soon. If you otherwise need help using this package, just let me know. I can be contacted at tyler.burns (at) drfz.de. Nonetheless, I think a visual tool like this is a great way to get some exposure to R programming (aside from enhancing your data analysis pipeline). 
