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


def methratio_parser(subparsers):
    workflow = subparsers.add_parser("wcg_methratio", 
        help = "Calculate each gene an average methylation ratio across all cells in promoter/genebody region")

    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--gene_bed", dest = "gene_bed", default = None,
        help = "Location of the reference genome bed file. ")
    group_input.add_argument("--cell_barcode", dest = "cell_barcode", default = None,
        help = "Location of the cell barcode list(generate by Preprocess snakemake pipeline). "
        "Cells which passed quality check.")
    group_input.add_argument("--peak_reference", dest = "peak_reference", default = None,
    	help = "Path to WCG.uniq.peak")
    group_input.add_argument("--meth_matrix", dest = "meth_matrix", default = None,
    	help = "Path to WCG.site_peak.h5")
 
    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("--outdir", dest = "out_dir", default = ".", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: current directory. ")
    group_output.add_argument("--outprefix", dest = "out_prefix", default = "MultiSpace", 
        help = "Prefix of output files. DEFAULT: MultiSpace. ")

    group_part = workflow.add_argument_group("Part arguments")
    group_part.add_argument("--region", dest = "region", choices = ['promoter', 'genebody'],default = "promoter", 
        help = "Type of methylation region. promoter or genebody.  "
        "If not specified, MultiSpace will use promoter as default. ")
    group_part.add_argument("--distance", dest = "distance", type = int, default = 2000,
        help = "Distance of gene promoter region. GENEBODY NOT REQUIRED! For example, 10000. "
        "If not specified, MultiSpace will take 2000 as default. ")
 


def WCGMethyMatrix(peak_reference, gene_bed, meth_matrix, cell_barcode, out_dir, out_prefix, region, distance):

	WCG_site = pd.read_table(peak_reference, header = None, names = ["chr","start","end"], delimiter = "_")	
	peak_list = open(peak_reference, 'rt')
	WCG_mat = h5py.File(meth_matrix,'r')
	# mat = WCG_mat['Mcsc']
	WCG_site_matrix = sparse.csc_matrix((WCG_mat['Mcsc']['data'][:],WCG_mat['Mcsc']['indices'][:],WCG_mat['Mcsc']['indptr'][:]), WCG_mat['Mcsc'].attrs['shape'])
	# WCG_site_matrix.index = WCG_site[0]
	# WCG_site_matrix.columns = barcode_list[0]


	if region == "genebody":
		print("Calculating" + " "+ str(region) + " "+ "methylation ratio...")
		WCG_mat = WCGGenebodyScore(WCG_site_matrix, gene_bed, peak_list, cell_barcode)
		WCG_mat.to_csv(os.path.join(out_dir + out_prefix + "_WCG_" +region + "_ratio.csv"),index = True)
	elif region == "promoter":
		distance = distance
		# print(distance)
		print("Calculating" + " "+ str(region) + " "+ "methylation ratio...")
		WCG_mat = WCGPromoterScore(distance, WCG_site_matrix, gene_bed, peak_list, cell_barcode)
		WCG_mat.to_csv(os.path.join(out_dir + out_prefix + "_WCG_" + str(distance) + "_" +region + "_ratio.csv"),index = True)



def ExtractGeneInfo(gene_bed):
	"""Extract gene information from gene bed file."""

	bed = pd.read_csv(gene_bed, sep="\t", header=None,index_col = False, names = ["name","transcript_id","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","Score","name2","orfstart","orfend","orfcodon_frame"])
	bed['tss'] = bed.apply(lambda x: x['txStart'] if x['strand']=='+' else x['txEnd'], axis=1)#tss on both -/+ strand，promoter = tss+2000/tss-2000

	bed['promoter'] = bed.apply(lambda x: tuple([x['tss']-2000, x['tss']+2000]), axis=1)
	bed['uid'] = bed.apply(lambda x: "%s@%s@%s"%(x['name2'], x['txStart'], x['txEnd']), axis=1)
	bed = bed.drop_duplicates(subset='uid', keep="first")
	gene_info = []
	for irow, x in bed.iterrows():
		gene_info.append([x['chrom'], x['txStart'], x['txEnd'], x['tss'],1, x['uid']])
	### [chrom_0, txstart_1, txend_2, tss_3,name_4, 1_5, uid_6]
	return(gene_info)


def AnnoPeak(gene_bed, peak_list):
	"""Generate a sparse matrix to store each cell peak site """

	peaks_info = []
	genes_info = ExtractGeneInfo(gene_bed)
	genes_list = []
	genes_info_tss = list()
	genes_info_full = list()


	for igene in range(len(genes_info)):
		tmp_gene = genes_info[igene]
		genes_list.append(tmp_gene[-1])
		genes_info_full.append(tmp_gene + [igene]) #[chrom_0, txstart_1, txend_2, tss_3, 1_4, uid_5, [igene]]
	### add index at the end of gene symbol
	genes = list(set([i.split("@")[0] for i in genes_list]))

	for ipeak, peak in enumerate(peak_list):
		peaks_tmp = peak.strip().split("_")
		peak = peak.replace("\n","")
		peaks_info.append([peaks_tmp[0], int(peaks_tmp[1]), int(peaks_tmp[1]), int(peaks_tmp[2]), 0, peak, ipeak])
		# peaks_info [chrom_0, center_1, start_2, end_3, 0_4, uid_5, [ipeak]]

	return(genes_info_full, peaks_info, genes_list, genes)



def WCGPromoterScore(distance, WCG_site_matrix, gene_bed, peak_list, cell_barcode):
	"""Calculate promoter region an average methylation level in each gene in each cell """

	barcode_list = pd.read_csv(cell_barcode,header = None)
	(genes_info_full, peaks_info, genes_list, genes) = AnnoPeak(gene_bed, peak_list)
	gene_distance = distance
	genes_peaks_score_array = sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)

	w = genes_info_full + peaks_info
	A = {}

	w.sort()
	for elem in w: 
		if elem[4] == 1: 
			A[elem[-1]] = [elem[0], elem[3]] #if gene, A[igene] = chr, tss
		else: #judge peak position
			dlist = []
			for gene_name in list(A.keys()): 
				g = A[gene_name] # g = chr start (genes)
				tmp_distance = elem[1] - g[1] # tmp_dis = peak site - gene tss
				if (g[0] != elem[0]) or tmp_distance > gene_distance: # compare chr(g[0]) and enhancer dis
					dlist.append(gene_name) 
				else:
					genes_peaks_score_array[gene_name, elem[-1]] = 1
			for gene_name in dlist:
				del A[gene_name]

	w.reverse()
	for elem in w: 
		if elem[4] == 1: 
			A[elem[-1]] = [elem[0], elem[3]]
		else: #judge peak position
			dlist = []
			for gene_name in list(A.keys()): 
				g = A[gene_name] # g = chr start (genes)
				tmp_distance = g[1] - elem[1] # tmp_dis = peak site - gene tss
				if (g[0] != elem[0]) or tmp_distance > gene_distance: # compare chr(g[0]) and enhancer dis
					dlist.append(gene_name) 
				else:
					genes_peaks_score_array[gene_name, elem[-1]] = 1
			for gene_name in dlist:
				del A[gene_name]


	genes_peaks_score_array_csr = sparse.csr_matrix(genes_peaks_score_array)

	#mean methy level of each cell in each gene
	row_sums = np.array(sparse.csr_matrix(genes_peaks_score_array_csr).sum(axis = 1))[:,0] 
	row_indices, col_indices = genes_peaks_score_array_csr.nonzero()
	genes_peaks_score_array_csr.data /= row_sums[row_indices]

	genes_peaks_score_array = genes_peaks_score_array_csr.tocsc()
	genes_cells_score_csc = genes_peaks_score_array.dot(WCG_site_matrix)

	score_cells_dict = {}
	score_cells_sum_dict = {}

	for igene, gene in enumerate(genes_list):
		score_cells_dict[gene] = igene
		score_cells_sum_dict[gene] = genes_cells_score_csc[igene, :].sum()


	score_cells_dict_dedup = {}
	score_cells_dict_max = {}
	for gene in genes:
		score_cells_dict_max[gene] = float("-inf")


	#remove duplicate genes
	for gene in genes_list:
		symbol = gene.split("@")[0]
		if score_cells_sum_dict[gene] > score_cells_dict_max[symbol]:
			score_cells_dict_dedup[symbol] = score_cells_dict[gene]
			score_cells_dict_max[symbol] = score_cells_sum_dict[gene]
	
	gene_symbol = sorted(score_cells_dict_dedup.keys())
	matrix_row = []
	for gene in gene_symbol:
		matrix_row.append(score_cells_dict_dedup[gene])
	score_cells_matrix = genes_cells_score_csc[matrix_row, :]

	WCG_promoter = pd.DataFrame(score_cells_matrix.toarray(),columns = barcode_list.loc[:,0].tolist(),index = genes)
	return(WCG_promoter)



def WCGGenebodyScore(WCG_site_matrix, gene_bed, peak_list,cell_barcode):
	"""Calculate genebody region an average methylation level in each gene in each cell """

	barcode_list = pd.read_csv(cell_barcode,header = None)
	(genes_info_full, peaks_info, genes_list, genes) = AnnoPeak(gene_bed, peak_list)
	genes_peaks_score_array = sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)

	w = genes_info_full + peaks_info
	A = {}

	w.sort()
	for elem in w: 
		if elem[4] == 1: 
			A[elem[-1]] = [elem[0], elem[1], elem[2]] #if gene, A[igene] = chr, txStart, txEnd
		else: #judge peak position
			dlist = []
			for gene_name in list(A.keys()): 
				g = A[gene_name] # g = chr start end(genes)
				if (g[0] == elem[0]) and (elem[1]>=g[1]) and (elem[1]<=g[2]) : # compare chr(g[0]) and peak site & gene txstart txend
					genes_peaks_score_array[gene_name, elem[-1]] = 1
				else:
					dlist.append(gene_name)
			for gene_name in dlist:
				del A[gene_name]
	#[chrom_0, center_1, start_2, end_3, 0_4, uid_5, [ipeak]]


	genes_peaks_score_array_csr = sparse.csr_matrix(genes_peaks_score_array)

	row_sums = np.array(sparse.csr_matrix(genes_peaks_score_array_csr).sum(axis = 1))[:,0]
	row_indices, col_indices = genes_peaks_score_array_csr.nonzero()
	genes_peaks_score_array_csr.data /= row_sums[row_indices]

	genes_peaks_score_array = genes_peaks_score_array_csr.tocsc()
	genes_cells_score_csc = genes_peaks_score_array.dot(WCG_site_matrix)

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
			score_cells_dict_dedup[symbol] = score_cells_dict[gene]
			score_cells_dict_max[symbol] = score_cells_sum_dict[gene]

	gene_symbol = sorted(score_cells_dict_dedup.keys())
	matrix_row = []
	for gene in gene_symbol:
		matrix_row.append(score_cells_dict_dedup[gene])
	score_cells_matrix = genes_cells_score_csc[matrix_row, :]

	WCG_genebody = pd.DataFrame(score_cells_matrix.toarray(),columns = barcode_list.loc[:,0].tolist(),index = genes)
	return(WCG_genebody)




   




