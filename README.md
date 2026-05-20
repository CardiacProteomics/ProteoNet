# ProteoNet

<!-- badges: start -->

<!-- Add badges here -->

<!-- badges: end -->

`ProteoNet` is an R package for network-based analysis and visualization of proteomics data. The package streamlines the construction, annotation, and visualization of protein interaction networks from differential proteomics experiments, with built-in support for clustering, overrepresentation analysis, and publication-ready plotting.

The package was developed for proteomics-focused systems biology workflows, with a particular emphasis on longitudinal and cardiac proteomics applications.

---

# Features

* Construction of protein-protein interaction networks from proteomics datasets
* Integration with STRING interaction data
* Automatic filtering and thresholding of interaction networks
* Detection of network subclusters using graph-based clustering algorithms
* Overrepresentation analysis (ORA) using user specified databases (e.g. Gene Ontology)
* Filters out redundancy of significant pathway terms
* Publication-ready network visualizations using `ggraph`
* Support for singleton protein integration and network augmentation

---

# Installation

## Install from GitHub

```r
# Install remotes if needed
install.packages("remotes")

# Install ProteoNet
remotes::install_github("CardiacProteomics/ProteoNet")
```

## Load the package

```r
library(ProteoNet)
```

---

# Workflow overview

A typical `ProteoNet` workflow consists of:

1. Preparing a list of proteins with associated statistics
2. Constructing a protein interaction network
3. Filtering interactions by confidence score
4. Detecting network clusters
5. Performing pathway enrichment analysis
6. Removing redundant pathways
7. Visualizing the resulting network

---

# Example workflow

```r
# Run the full ProteoNet analysis pipeline
res <- proteonet_pipeline(
  reference,          # Label for saving output
  my_geneset,         # Vector of genes/proteins to include in the network
  species,            # Species identifier (9606 = human, 10090 = mouse)
  min_cluster_size,   # Minimum number of proteins required for a subcluster
  score_threshold,    # STRING interaction confidence threshold
  selection,          # Can only be "fdr" (for picking representative out of redundant pathways)
  databases_tested,   # Databases used for overrepresentation analysis (e.g. "gocc_human", "gobp_human", "gocc_mouse")
  ora_min,            # Minimum pathway size for enrichment analysis
  ora_max,            # Maximum pathway size for enrichment analysis
  folder_string,      # Main output directory
  folder_genesets,    # Directory for exported gene set files
  folder_results,     # Directory for result tables
  folder_figures,     # Directory for generated figures
  universe,           # Background gene set for enrichment analysis
  threshold_mean,     # Mean pathway cluster correlation threshold 
  threshold_min,      # Minimum pathway cluster correlation threshold 
  min_score = 0.1       # Minimum node connectivity score
)

# Replot the generated network with custom figure settings
res2 <- replot_network(
  res$network,                    # Network object returned by proteonet_pipeline()
  folder_figures,                # Folder where figures should be saved
  network_height = 20,           # Height of exported figure
  network_width = 20,            # Width of exported figure
  label_size = 3,                # Size of node labels
  legend_size = 4,               # Size of legend text
  unlist(res$labels),            # Labels extracted from pipeline output
  reference                      # Original reference data frame
)
```



---

# Network construction

Networks are generated using known or predicted protein-protein interactions. Interactions can be filtered using confidence score thresholds to control network density.

Clustering is performed using graph-based community detection methods to identify biologically meaningful subnetworks.

---

# Pathway enrichment analysis

`ProteoNet` includes support for overrepresentation analysis (ORA) using Gene Ontology databases.

Supported annotation resources may include:

* GO Biological Process
* GO Cellular Component
* Additional enrichment databases depending on configuration

The most significantly enriched pathways can be automatically annotated directly onto network visualizations.

---

# Visualization

The package uses `ggplot2` and `ggraph` for customizable and publication-quality network visualizations.

Visualization features include:

* Cluster-aware network layouts
* Node coloring by cluster
* Highlighting of singleton proteins
* Edge styling based on interaction type
* Automatic pathway labeling

---

# Output

The pipeline returns structured output objects containing:

* Network graphs
* Cluster assignments
* Enrichment analysis results
* Layout coordinates
* Visualization-ready data frames

This makes it easy to perform downstream analyses or generate custom figures.

---

# Citation

If you use `ProteoNet` in your work, please cite:

```text
Citation information coming soon.
```

---

# Development status

`ProteoNet` is under active development. Contributions, feature requests, and bug reports are welcome.


# Issues and bug reports

If you encounter a bug or have a feature request, please open an issue on GitHub.

---

# License

Specify your package license here.

For example:

```text
MIT License
```
