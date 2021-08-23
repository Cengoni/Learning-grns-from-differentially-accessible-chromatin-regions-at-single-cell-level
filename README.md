This repository contains large part of the scripts and notebooks used during my master thesis project. Most include the analyses described in the report and minor utility scripts for visualization and processing.
In the Pre-processing folder there are:
- preprocessing.R, which includes the pre-processing of chromatin accessibility data, largely conducted using Seurat R package. These computational steps correspond to the ones described in 3.1.2 Quality control of cells and fragments and replicates integration, 3.2.1 Clustering and 3.2.2 Cell identity annotation. The input to this file are the pre-processed peaks-by-cell matrices from CellRanger. 
- peaks_comparison.R, which includes the analyses reported in 3.1.1 Peak calling: comparison between MACS2 and CellRanger.


What is NOT included are the earlier steps performed in CellRanger and MACS2, which are conducted as reported in the Material and Methods 5.1.1 Reads pre-processing and 5.1.2 Peak calling and matrix formation from command line.

In the Downstream_analysis folder there are:
- network_analyses.ipynb, which includes network filtering by gene-body openness and by expression data, with relative statistical tests and comparison between scATAC and snucRNA. These correspond to the analyses described in 3.3.2 Gene-body openness as a filtering criterion for TFs and 3.3.4 Comparison with snucRNA-seq at TFs and TGs level.
- celltype_specificity.ipynb, which reports the comparison between TF and TGs between the three cell-types. This includes the plots used to create Figure 3.8.
- functional_analyses.ipynb, which includes all the plots and analyses performed using functional annotation and enrichment of GO terms. These include the analyses discussed in chapter 3.4 Functional analysis and 3.2.3 Functional enrichment of clusters confirms their biological annotation. 
- other_figures.R, which includes the code for figure 3.5. C and figure 3.9 A.
- extract_tfs_position.R, which is an utility script used to extract TFs coordinates (1kb upstream), required for gene-body openness filtering (network_analyses.ipynb)

What is NOT included is the network inference step, which includes motif collection and motif mapping, conducted using the tools mentioned in 5.3.1 Motif collection and motif mapping, and the MINI-AC network inference pipeline, as it is an unpublished method of VIB Comparative Network Biology (CNB) group of Prof. DR. Klaas Vandepoele.
