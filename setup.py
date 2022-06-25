#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   setup.py
@Time    :   2021/04/16 12:34:01
@Author  :   Zhanhe Chang
@Contact :   changzhanhe98@163.com
@License :   (C)Copyright 2021-2022, Zhanhe Chang
'''


import sys,os

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

exec(open('MultiSpace/version.py').read())

with open("requirements.txt") as f:
    requirements = f.read().splitlines()



def main():
    setup(
        name = "MultiSpace",
        version = __version__,
        # package_dir = {'':'MultiSpace'},
        packages = ['MultiSpace', 'MultiSpace.Generatescorematrix', 'MultiSpace.MappingscCell', 'MultiSpace.MultiSpace_init', ],
        package_data = {'MultiSpace':['Snakemake/modules/*', 'Snakemake/modules/scripts/*', 'Snakemake/*', 'annotations/*']},
        install_requires = requirements,
        setup_requires = requirements,
        entry_points={
            'console_scripts': [
            'MultiSpace=MultiSpace.start:main'
            ]
        },
        include_package_data = True,
        author = "Zhanhe Chang",
        author_email = "changzhanhe98@163.com",
        description = "MultiSpace (Single-cell Multi Omics Analysis In Space) is a pipeline which can map single cell to spatial location and get spatially epigenetic signal. ",
        license = "GPL-3.0",
        url = "https://github.com/Changzhanhe/MultiSpace",
        classifiers = [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        python_requires=">=3.7",
    )

if __name__ == "__main__":
    main()


