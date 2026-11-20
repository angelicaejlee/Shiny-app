# app.R

**app.R** is a k-mean shiny app to run and visualize k-means clustering by PCA, UMAP, and NMF plots.

On the left side panel, there is a menu to choose the method including PCA, UMAP, and NMF.
The user can also choose the number of clusters from minimum of 1 cluster to maximum of 9 clusters.
The main panel includes the plot and a table that shows clusters vs. TB status and adjusted rand index.
Different clusters are distinguished by different colors and the centers of each cluster are marked by x.
The shapes of each data points differ depending on the TB status.
R program is used to run this app.
