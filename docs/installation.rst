.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------




System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* Python >= 3.7


We recommend to create an independent conda environment for MultiSpace. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Install the stable version
>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for MultiSpace.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   git clone https://github.com/Changzhanhe/MultiSpace.git
   cd MultiSpace
   conda env create -f environment.yml -n multispace
   conda activate multispace

Step 2 Install MultiSpace.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   python setup.py install