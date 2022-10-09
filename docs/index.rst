.. MultiSpace documentation master file, created by
   sphinx-quickstart on Sun Jun 19 10:33:38 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MultiSpace's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


|Docs| |Conda|

.. |Docs| image:: https://readthedocs.org/projects/multispace/badge/?version=latest
   :target: https://multispace.readthedocs.io

.. |Conda| image:: https://anaconda.org/changzhanhe/multispace/badges/version.svg


MultiSpace (Single-cell Multi-omic Analysis In Space), a computational framework that combines single-cell multi-omic data such as scCOOL-seq with spatial transcriptomic information. MultiSpace first projects single cells with multi-omic modalities into spatial locations by deconvolution-based cell mapping with transcriptomic data, thereby reconstructing RNA expression, DNA methylation, and chromatin accessibility information in space with high efficiency. MultiSpace then provided several pre-processing functions including generating gene/transposons expression and DNA methylation level, and calculating gene activity scores for both single-cell and spatial data. Besides, MultiSpace also provides rich downstream analysis functions including 1) single and multi-modality clustering, visualization, differential element calling, and trajectory detection, 2) spatial domain identification, spatial variable element calling, and visualization. All these components have been packed using Snakemake workflow for convenient distribution. 


.. image:: _static/img/workflow.png
   :width: 100%
   :align: center


.. include:: release_notes/0.0.1.rst


.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   usage
   examples
   release_notes/index


