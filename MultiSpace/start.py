#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse as ap

from MultiSpace.version import __version__
from MultiSpace.MultiSpace_init.MultiSpace_init import pipelineinit_parser,PipelineConfig
from MultiSpace.GeneMethyscore.GeneMethyscore import methratio_parser,WCGMethyMatrix
from MultiSpace.GeneActivity.GeneActivity import geneactivity_parser, GCHScoreMatrix
from MultiSpace.MappingscCell.MappingscCell import gettingepisignal_parser, Mapping


def main():
    """
    Add main function argument parsers.
    """

    parser = ap.ArgumentParser(description = "MultiSpace(Single-cell Multi Omics Analysis In Space) is a multi-omics pipeline integrated RNA Expression, DNA methylation and Chromatin Accessibility analysis built using snakemake.")
    parser.add_argument("-v", "--version", action = "store_true", help = "Print version info.")

    subparsers = parser.add_subparsers(dest = "subcommand")

    pipelineinit_parser(subparsers)
    methratio_parser(subparsers)
    geneactivity_parser(subparsers)
    gettingepisignal_parser(subparsers)

    logging.basicConfig(format="%(levelname)s: %(message)s", stream=sys.stderr)
    args = parser.parse_args()

    version = __version__

    if args.version:
        print(version)
        exit(0)
    elif args.subcommand == "pipeline_init":
        PipelineConfig(species = args.species, samplesheet = args.samplesheet, star_annotation = args.star_annotation, 
            fasta = args.fasta, lambda_fasta = args.lambda_fasta, fasta_fai = args.fasta_fai, star_index = args.star_index, directory = args.directory)
    elif args.subcommand == "gch_geneactivity":
        GCHScoreMatrix(file_path = args.file_path , gene_bed = args.gene_bed, cell_barcode = args.cell_barcode, out_dir = args.out_dir, 
            out_prefix = args.out_prefix, distance = args.distance)
    elif args.subcommand == "wcg_methratio":
        WCGMethyMatrix(peak_reference = args.peak_reference, gene_bed = args.gene_bed, meth_matrix = args.meth_matrix, cell_barcode = args.cell_barcode, 
            out_dir = args.out_dir, out_prefix = args.out_prefix, region = args.region, distance = args.distance)
    elif args.subcommand == "getting_episignal":
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



 



