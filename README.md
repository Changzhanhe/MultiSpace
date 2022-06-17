# MultiSpace
MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, DNA methylation and Chromatin Accessibility analysis built using snakemake.MultiSpace support scCOOL-seq, scNMT-seq for DNA analysis; Smart-seq2 for RNA analysis. MultiSpace combines several tools and packages to create an integrative pipeline, which enables three omics anaylsis from raw sequencing data (fastq files as input) through alignment, quality control, cell filtering, methylation site calling and filtering, generate site by cell and bin by cell matrix in single cell level. Besides preprocessing, MultiSpace also provides several downstream analysis functions, including (1) gene genebody/promoter DNA methylation ratio, (2) using the regulatory potential model to calculate gene activity score, (3) mapping single cell to spatial location and get spatial epigenetic (DNA methylation/Chromatin Accessibility) signal.

![avatar](docs/imag.png)

