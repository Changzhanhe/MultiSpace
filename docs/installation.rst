.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------



System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* Python >= 3.8


We recommend to create an independent conda environment for MultiSpace. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Install MultiSpace using conda(stable version)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for MultiSpace.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   # Create a python3.8 environment for installing MultiSpace.
   conda create -n multispace python=3.8
   # Install through the following commands:
   conda config --add channels defaults
   conda config --add channels dongqingsun
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --add channels changzhanhe

Step 2 Install MultiSpace.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   # To make the installation faster, we recommend using mamba
   conda install mamba -c conda-forge
   mamba install -c changzhanhe multispace


Install MultiSpace from Github(developing version)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

::::::::::::::::::::::::::::::::::::::::::::
:: 

   git clone https://github.com/Changzhanhe/MultiSpace.git
   cd MultiSpace
   # Create a conda environment for MultiSpace.
   conda env create -f environment.yml -n multispace
   # Installing package
   conda activate multispace
   python setup.py install






