# -*- coding: utf-8 -*-
# @Author: Zhanhe Chang

import os,sys
import tables
import h5py
import argparse

import numpy as np
import argparse as ap
import pandas as pd
import scipy.sparse as sparse

from pkg_resources import resource_filename
from pandas import Series,DataFrame


def geneactivity_parser(subparsers):
	workflow = subparsers.add_parser("gch_geneactivity", 
		help = "Calculate each gene an activity score across all cells")

	group_input = workflow.add_argument_group("Input arguments")
	group_input.add_argument("--gene_bed", dest = "gene_bed", default = None,
		help = "Location of the reference genome bed file. ")
	group_input.add_argument("--cell_barcode", dest = "cell_barcode", default = None,
		help = "Location of the cell barcode list(generate by Preprocess snakemake pipeline). "
		"Cells which passed quality check.")
	group_input.add_argument("--file_path", dest = "file_path", default = None,
		help = "Path to unipeak file and site_peak.h5 file")
 
	group_output = workflow.add_argument_group("Output arguments")
	group_output.add_argument("--out_dir", dest = "out_dir", default = ".", 
		help = "Path to the directory where the result file shall be stored. DEFAULT: current directory. ")
	group_output.add_argument("--out_prefix", dest = "out_prefix", default = "MultiSpace", 
		help = "Prefix of output files. DEFAULT: MultiSpace. ")

	group_part = workflow.add_argument_group("Part arguments")
	group_part.add_argument("--distance", dest = "distance", type = int, default = 10000,
		help = " Gene score decay distance, could be optional from 1kb (promoter-based regulation) to 10kb (enhancer-based regulation). "
		"DEFAULT: 10000.")
 


def get_peakpattern(pathdir):
	allfiles = os.listdir(pathdir)
	peakpattern = []
	pattern =  "^GCH.*uniq.peak$"
	for file in allfiles:
		if re.fullmatch(pattern, file):
			peakpattern.append(file.split(".")[1])
		else:
			pass
	return(peakpattern)



def GCHScoreMatrix(file_path , gene_bed, cell_barcode, out_dir, out_prefix, distance):

  
	chrpat = get_peakpattern(file_path)
	dflist = []
	barcode_list = pd.read_csv(cell_barcode,header = None)
	
	for chr,ichr in enumerate(chrpat):
		print("Calculating" + " "+ str(chrpat[chr]) + " "+ "activity score...")
		peak_reference = file_path + "GCH." + chrpat[chr] + ".uniq.peak"
		GCH_site = pd.read_table(peak_reference, header = None, names = ["chr","start","end"], delimiter = "_")	
		peak_list = open(peak_reference, 'rt')
		gch_matrix = file_path + "GCH." + chrpat[chr] + ".site_peak.h5"
		GCH_mat = h5py.File(gch_matrix,'r')
		GCH_site_matrix = sparse.csc_matrix((GCH_mat['Mcsc']['data'][:],GCH_mat['Mcsc']['indices'][:],GCH_mat['Mcsc']['indptr'][:]), GCH_mat['Mcsc'].attrs['shape'])


		chrlist = []
		if str(ichr).split("_")[0] == "chr1":
			for i in range(1,7):
				chrlist.append("chr" + str(i))
		elif str(ichr).split("_")[0] == "chr7":
			for i in range(7,13):
				chrlist.append("chr" + str(i))
		elif str(ichr).split("_")[0] == "chr13":
			for i in range(13,20):
				chrlist.append("chr" + str(i))
			chrlist.extend(['chrM','chrX','chrY'])

		
		score_df = "GCH_mat.%s" % ichr
		(genes_info, peaks_info, genes_list, genes) = AnnoPeak(gene_bed, peak_list, chrlist)
		score_df = GCHActivityScore(distance, GCH_site_matrix, genes_info, peaks_info, genes_list, genes, barcode_list)
		dflist.append(score_df)

		
	GCH_mat = pd.concat([x for x in dflist]).groupby(['Gene']).sum()
	GCH_mat.to_csv(os.path.join(out_dir + out_prefix + "_GCH_" + "score.csv"),index = True)

	

def AnnoPeak(gene_bed, peak_list, chrlist):
	"""Generate a sparse matrix to store each cell peak site """
	
	peaks_info = []
	genes_info = []
	genes_list = []


	fhd = open(gene_bed,'r')
	fhd.readline() 
	for line in fhd:
		line = line.strip().split('\t')
		if line[2] in chrlist:
			if not line[0].startswith('#'):
				if line[3] == "+":
					genes_info.append((line[2], int(line[4]), 1, "%s@%s@%s" % (line[12], line[2], line[4])))
				else:
					genes_info.append((line[2], int(line[5]), 1, "%s@%s@%s" % (line[12], line[2], line[5])))
					# gene_info [chrom, tss, 1, gene_unique]
	fhd.close()
	genes_info = list(set(genes_info))

	for igene in range(len(genes_info)):
		tmp_gene = list(genes_info[igene])
		genes_list.append(tmp_gene[3])
		tmp_gene[3] = igene #no.XX
		genes_info[igene] = tmp_gene #['chr1',171899539,1,32980]
	genes = list(set([i.split("@")[0] for i in genes_list]))

	
	for ipeak, peak in enumerate(peak_list):
		peaks_tmp = peak.strip().split('_')
		peaks_info.append([peaks_tmp[0], int(peaks_tmp[1]), 0, ipeak])
		
		
	return(genes_info, peaks_info, genes_list, genes)


def GCHActivityScore(distance, GCH_site_matrix, genes_info, peaks_info, genes_list, genes, barcode_list):

	gene_distance = 10000
	Sg = lambda x: 2**(-x)
	decay = float(gene_distance)
	gene_distance = 15 * decay

	genes_peaks_score_array = sparse.dok_matrix((len(genes_info), len(peaks_info)), dtype=np.float16)

	w = genes_info + peaks_info
	A = {}

	w.sort()
	for elem in w: 
		if elem[2] == 1: 
			A[elem[-1]] = [elem[0], elem[1]]
		else: 
			dlist = []
			for gene_name in list(A.keys()): 
				g = A[gene_name] # g = chr start (genes)
				tmp_distance = elem[1] - g[1] # tmp_dis = peak site - gene tss
				if (g[0] != elem[0]) or (tmp_distance > gene_distance): # compare chr and enhancer dis
					dlist.append(gene_name) 
				else:
					genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
			for gene_name in dlist:
				del A[gene_name]

	w.reverse()
	for elem in w:
		if elem[2] == 1:
			A[elem[-1]] = [elem[0], elem[1]]
		else:
			dlist = []
			for gene_name in list(A.keys()):
				g = A[gene_name]
				tmp_distance = g[1] - elem[1]
				if (g[0] != elem[0]) or tmp_distance > gene_distance:
					dlist.append(gene_name)
				else:
					genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
			for gene_name in dlist:
				del A[gene_name]

	genes_peaks_score_csc = genes_peaks_score_array.tocsc()
	genes_cells_score_csc = genes_peaks_score_csc.dot(GCH_site_matrix)

	score_cells_dict = {}
	score_cells_sum_dict = {}


	for igene, gene in enumerate(genes_list):
		score_cells_dict[gene] = igene
		score_cells_sum_dict[gene] = genes_cells_score_csc[igene, :].sum()


	score_cells_dict_dedup = {}
	score_cells_dict_max = {}

	for gene in genes:
		score_cells_dict_max[gene] = float("-inf")


	for gene in genes_list:
		symbol = gene.split("@")[0]
		if score_cells_sum_dict[gene] > score_cells_dict_max[symbol]:
			score_cells_dict_dedup[symbol] = score_cells_dict[gene] #in 1500000 range consider all peak
			score_cells_dict_max[symbol] = score_cells_sum_dict[gene]
	gene_symbol = sorted(score_cells_dict_dedup.keys())
	matrix_row = []
	for gene in gene_symbol:
		matrix_row.append(score_cells_dict_dedup[gene])

	score_cells_matrix = genes_cells_score_csc[matrix_row, :]

	GCH_score = pd.DataFrame(score_cells_matrix.toarray(),columns = barcode_list.loc[:,0].tolist(),index = genes)
	GCH_score.index.name = 'Gene'
	return(GCH_score)





























