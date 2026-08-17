# MSc_IRP_229004366
## Software environment and reproducibility

Bulk RNA-seq preprocessing was performed on the University of Leicester ALICE high-performance computing cluster using STAR v2.7.11a, SAMtools v1.17 and LiBiNorm v2.5. Reads were aligned to the GRCh38 primary assembly using an Ensembl release 116 annotation (`Homo_sapiens.GRCh38.116.gtf`).

Downstream analysis was performed locally using R v4.6.1 on Linux Mint 22.1 through RStudio v2026.07.0. Key R packages included DESeq2 v1.52.0, dplyr v1.2.1, ggplot2 v4.0.3, doRothEA v1.23.0, ggVennDiagram v1.5.7, ggrepel v0.9.8, igraph v2.3.3 and ggraph v2.2.2. Full R session information, including all loaded package versions and library paths, is provided in `session_info.txt`.
