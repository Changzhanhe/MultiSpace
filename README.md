# MultiSpace

![Documentation Status](https://readthedocs.org/projects/multispace/badge/?version=latest)
[![Anaconda-Server Badge](https://anaconda.org/changzhanhe/multispace/badges/license.svg)](https://anaconda.org/changzhanhe/multispace)
[![Anaconda-Server Badge](https://anaconda.org/changzhanhe/multispace/badges/installer/conda.svg)](https://conda.anaconda.org/changzhanhe)
[![Anaconda-Server Badge](https://anaconda.org/changzhanhe/multispace/badges/platforms.svg)](https://anaconda.org/changzhanhe/multispace)


MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, DNA methylation and Chromatin Accessibility analysis built using snakemake.MultiSpace support scCOOL-seq, scNMT-seq for DNA analysis; Smart-seq2 for RNA analysis. MultiSpace combines several tools and packages to create an integrative pipeline, which enables three omics anaylsis from raw sequencing data (fastq files as input) through alignment, quality control, cell filtering, methylation site calling and filtering, generate site by cell and bin by cell matrix in single cell level. Besides preprocessing, MultiSpace also provides several downstream analysis functions, including (1) gene genebody/promoter DNA methylation ratio, (2) using the regulatory potential model to calculate gene activity score, (3) mapping single cell to spatial location and get spatial epigenetic (DNA methylation/Chromatin Accessibility) signal.

![avatar](docs/_static/img/workflow.png)

## Documentation
For full installation and usage of MultiSpace, please refer to the [documentation](https://multispace.readthedocs.io/en/latest/).


## Change Log
### v0.0.1
* Release MultiSpace.
* Use Snakemake to preprocess raw data. Add pipeline initiation function.
* Add mapping single cell to spatial location, DNA methylation ratio, gene activity score using RP-model functions.

## Installation
### Use the following commands to install Minicoda3：
``` bash
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
### Install MultiSpace
```bash
git clone https://github.com/Changzhanhe/MultiSpace.git
cd MultiSpace
conda env create -f environment.yml -n multispace
conda activate multispace
python setup.py install
```

## Usage
```bash
usage: MultiSpace [-h] [-v] {Pipelineinit,Scorematrix,Mappingcell} ...

MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, DNA methylation and Chromatin Accessibility analysis built using snakemake.

positional arguments:
  {Pipelineinit,Scorematrix,Mappingcell}
    Pipelineinit        Initialize the MultiSpace preprocessing workflow in a given directory. This will install the snakemake rules and a config file in this directory. You can configure the config file
                        according to your needs, and run the workflow with Snakemake.
    Scorematrix         Calculate each gene a gene by cell score matrix across all cells. WCG: Genebody/Promoter methylation ratio matrix. GCH: Gene activity score matrix.
    Mappingcell         Map single cell to spatial location according to expression similarity and get spatial epigenetic signal.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version info.
```
