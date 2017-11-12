## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE, warning = FALSE)
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

## ------------------------------------------------------------------------
library(Sconify)
library(ggplot2)
set.seed(12043)
final <- post.processing(scone.output = scone.output,
                         cell.data = combined,
                         input = input.markers)
combined # input data
scone.output # scone-generated data
final # the data after post-processing

# tSNE map shows highly responsive population of interest
qplot(final$`bh-SNE1`, 
      final$`bh-SNE2`, 
      color = final$`pSTAT5(Nd150)Di.IL7.change`, 
      xlab = "bh-SNE1",
      ylab = "bh-SNE2") + 
    labs(color = "IL7 -> pSTAT5 change") + 
    scale_color_gradientn(colors = c("black", "yellow")) 

# tSNE map now colored by q value
qplot(final$`bh-SNE1`, 
      final$`bh-SNE2`, 
      color = final$`pSTAT5(Nd150)Di.IL7.qvalue`, 
      xlab = "bh-SNE1",
      ylab = "bh-SNE2") + 
    labs(color = "IL7 -> pSTAT5 -log10(qvalue)") + 
    scale_color_gradientn(colors = c("black", "yellow"))

