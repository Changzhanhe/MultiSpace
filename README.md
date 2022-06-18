# MultiSpace
MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, DNA methylation and Chromatin Accessibility analysis built using snakemake.MultiSpace support scCOOL-seq, scNMT-seq for DNA analysis; Smart-seq2 for RNA analysis. MultiSpace combines several tools and packages to create an integrative pipeline, which enables three omics anaylsis from raw sequencing data (fastq files as input) through alignment, quality control, cell filtering, methylation site calling and filtering, generate site by cell and bin by cell matrix in single cell level. Besides preprocessing, MultiSpace also provides several downstream analysis functions, including (1) gene genebody/promoter DNA methylation ratio, (2) using the regulatory potential model to calculate gene activity score, (3) mapping single cell to spatial location and get spatial epigenetic (DNA methylation/Chromatin Accessibility) signal.

![avatar](docs/imag.png)

## Change Log
### v0.0.1
* Release MultiSpace.
* Please prepare input files and folder structure following the installation instructions.
* Add mapping single cell to spatial location, DNA methylation level, chromatin accessibility level functions.

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
usage: MultiSpace [-h] [-v] {pipeline_init,wcg_methratio,gch_geneactivity,getting_episignal} ...

MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, 
DNA methylation and Chromatin Accessibility analysis built using snakemake.

positional arguments:
  {pipeline_init,wcg_methratio,gch_geneactivity,getting_episignal}
    pipeline_init       Initialize the MultiSpace preprocessing workflow in a given directory. 
                        This will install the snakemake rules and a config file in this directory. You can configure the config file
                        according to your needs, and run the workflow with Snakemake
    wcg_methratio       Calculate each gene an average methylation ratio across all cells in promoter/genebody region
    gch_geneactivity    Calculate each gene an activity score across all cells
    getting_episignal   Map single cell to spatial location and get spatial epigenetic signal.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version info.
```
