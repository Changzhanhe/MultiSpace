# -*- coding: utf-8 -*-
# @Author: Zhanhe Chang

import os
import shutil
import argparse as ap
from jinja2 import Template
from pkg_resources import resource_filename

def pipelineinit_parser(subparsers):
	"""
	Add main function snakemake preprocess pipeline argument parsers.
	"""
	workflow = subparsers.add_parser("Pipelineinit", help = "Initialize the MultiSpace preprocessing workflow in a given directory. "
		"This will install the snakemake rules and a config file in this directory. "
		"You can configure the config file according to your needs, and run the workflow with Snakemake ")

	# Input files arguments
	group_input = workflow.add_argument_group("Input files arguments")
	group_input.add_argument("--species", dest = "species", default = "mm10", 
		choices = ["hg38", "mm10"], 
		help = "Specify the genome assembly (hg38 for human and mm10 for mouse). DEFAULT: mm10.")
	group_input.add_argument("--samplesheet", dest = "samplesheet", default = "", type = str, 
		help = "Path to sample names stored in a sheet."
		"Row: sample name.")


	# Running files arguments
	group_runningdir = workflow.add_argument_group("Running and output files arguments")
	group_runningdir.add_argument("--directory", dest = "directory", default = "MultiSpace", type = str,
		help = "Path to the directory where the workflow shall be initialized and results shall be stored. DEFAULT: MultiSpace."
		"Path to where the config.yaml is stored.")


	# Reference genome arguments
	group_reference = workflow.add_argument_group("Reference genome arguments")
	group_reference.add_argument("--fasta", dest = "fasta",  type = str, 
		required = True,
		default = "",
		help = "Genome fasta file for mapping."
		"Users should provide fasta file for human and mouse only containing chrN, where N is the name of defined chromosome.")
	group_reference.add_argument("--fasta_fai", dest = "fasta_fai", type = str,
		required = True,
		default = "",
		help = "Genome fasta file index."
		"User can create fasta.fai using samtools faidx.")
	group_reference.add_argument("--lambda_fasta", dest = "lambda_fasta", type = str,
		required = True,
		default = "",
		help = "Genome fasta file containing lambda sequence for bsmap mapping."
		"Users can add lambda sequence to fasta file showed upper.")
	group_reference.add_argument("--rna_annotation", dest = "rna_annotation", type = str,
		required = True,
		default = "",
		help = "Path of the annotation file required for . "
		"Users can define or provide annotation files themselves.")
	group_reference.add_argument("--star_index", dest = "star_index", type = str,
		default = "",
		help = "Path of the reference index file for  STAR mapping."
		"Users need to build the index file for the reference using command "
		"STAR --runThreadN N --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ref.fa --sjdbGTFfile refGene.gtf")


def PipelineConfig(species, samplesheet, rna_annotation, fasta, fasta_fai, lambda_fasta, star_index, directory):
	"""
	Generate snakemake preprocess config.yaml file.
	"""

	try:
		os.makedirs(directory)
	except OSError:
		# either directory exists (then we can ignore) or it will fail in the
		# next step.
		pass

	pkgpath = resource_filename('MultiSpace', 'Snakemake')
	annopath = resource_filename('MultiSpace', 'annotations')
	template_file = os.path.join(pkgpath, "config_template.yaml")
	repeatLINEfile = os.path.join(annopath, species + ".repeat.LINE.merge.bed")
	repeatLTRfile = os.path.join(annopath, species + ".repeat.LTR.merge.bed")
	DNAdir = os.path.join(directory, "DNA/")
	RNAdir = os.path.join(directory, "RNA/")
	configfile = os.path.join(directory, "config.yaml")
	config_template = Template(open(template_file, "r").read(), trim_blocks=True, lstrip_blocks=True)
	with open(configfile, "w") as configout:
			configout.write(config_template.render(
				#running and output file arguments
				directory = os.path.abspath(directory),

				#input files arguments
				species = species,
				samplesheet = os.path.abspath(samplesheet),

				#reference genome arguments
				fasta = os.path.abspath(fasta),
				fasta_fai = os.path.abspath(fasta_fai),
				lambda_fasta = os.path.abspath(lambda_fasta),
				rna_annotation = os.path.abspath(rna_annotation),
				star_index = os.path.abspath(star_index),
				repeat_line = os.path.abspath(repeatLINEfile),
				repeat_ltr = os.path.abspath(repeatLTRfile)))
				# dnadir = os.path.abspath(DNAdir),
				# rnadir = os.path.abspath(RNAdir)))


	source = os.path.join(pkgpath, "Snakefile")
	modules = os.path.join(pkgpath, "modules")
	target = os.path.join(directory, "Snakefile")
	dest = os.path.join(directory, "modules")
	shutil.copy(source, target)
	shutil.copytree(modules, dest)


