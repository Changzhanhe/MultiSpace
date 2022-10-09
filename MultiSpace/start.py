#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse as ap

from MultiSpace.version import __version__
from MultiSpace.MultiSpace_init.MultiSpace_init import pipelineinit_parser,PipelineConfig
from MultiSpace.Generatescorematrix.Generatescorematrix import generatescore_parser,Generate_scorematrix
from MultiSpace.MappingscCell.MappingscCell import mappingcelltospatial_parser, Mapping


def main():
    """
    Add main function argument parsers.
    """

    parser = ap.ArgumentParser(description = "MultiSpace(Single-cell Multi-Omics Analysis In Space) is a s a computational framework that combines single-cell multi-omic data such as scCOOL-seq with spatial transcriptomic information.")
    parser.add_argument("-v", "--version", action = "store_true", help = "Print version info.")

    subparsers = parser.add_subparsers(dest = "subcommand")

    pipelineinit_parser(subparsers)
    generatescore_parser(subparsers)
    mappingcelltospatial_parser(subparsers)

    logging.basicConfig(format="%(levelname)s: %(message)s", stream=sys.stderr)
    args = parser.parse_args()

    version = __version__

    if args.version:
        print(version)
        exit(0)
    elif args.subcommand == "Pipelineinit":
        PipelineConfig(species = args.species, samplesheet = args.samplesheet, star_annotation = args.star_annotation, 
            fasta = args.fasta, lambda_fasta = args.lambda_fasta, fasta_fai = args.fasta_fai, star_index = args.star_index, directory = args.directory)
    elif args.subcommand == "Scorematrix":
        Generate_scorematrix(file_path = args.file_path , species = args.species, cell_barcode = args.cell_barcode, out_dir = args.out_dir, matrixtype = args.matrixtype,
            out_prefix = args.out_prefix, distance = args.distance, region = args.region)
    elif args.subcommand == "Mappingcell":
        Mapping(sc_count_file = args.sc_count_file , sc_anno_file = args.sc_anno_file, st_count_file = args.st_count_file, sc_scale_factor = args.sc_scale_factor, 
            st_scale_factor = args.st_scale_factor, out_dir = args.out_dir, out_prefix = args.out_prefix, 
            model_dir = args.model_dir, normalize = args.normalize, gene_use = args.gene_use, ntopics_list = args.ntopics_list, spatial_location = args.spatial_location,
            epi_binfile = args.epi_binfile, epi_feature = args.epi_feature)
    else:
        parser.print_help()
        exit(1)
    exit(0)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)



 



