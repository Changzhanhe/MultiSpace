# -*- coding: utf-8 -*-
# @Author: Zhanhe Chang
import os,re
import h5py
import scipy
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.spatial import distance_matrix
from collections import defaultdict

from STRIDE.Deconvolution import *


def mappingcelltospatial_parser(subparsers):
	workflow = subparsers.add_parser("Mappingcell", 
		help = "Map single cell to spatial location and get spatial epigenetic signal.")

	group_input = workflow.add_argument_group("Input arguments")
	group_input.add_argument("--sc_count_file", dest = "sc_count_file", default = None,
		help = "Location of the single-cell count matrix file. "
		"It could be '.h5' or tab-separated plain-text file with genes as rows and cells as columns. ")
	group_input.add_argument("--sc_celltype_file", dest = "sc_anno_file", default = None,
		help = "Location of the single-cell celltype annotation file. "
		"The file should be a tab-separated plain-text file without header. "
		"The first column should be the cell name, and the second column should be the corresponding celltype labels. ")
	group_input.add_argument("--st_count_file", dest = "st_count_file", default = None,
		help = "Location of the spatial gene count file. "
		"It could be '.h5' or tab-separated plain-text file with genes as rows and spots as columns. ")
	group_input.add_argument("--gene_use", dest = "gene_use", default = None,
		help = "Location of the gene list file used to train the model. "
		"It can also be specified as 'All', but it will take a longer time. "
		"If not specified, MultiSpace will find differential marker genes for each celltype, and use them to run the model. ")
	group_input.add_argument("--spatial_location", dest = "spatial_location", default = None,
		help = "Location of tissue spatial coordinates")
	group_input.add_argument("--model_dir", dest = "model_dir", default = None,
		help = "If users have the pre-trained model using the same scRNA-seq dataset, please provide the path of 'model' directory.")
	group_input.add_argument("--epi_binfile", dest = "epi_binfile", default = None,
		help = "Location of WCG/GCH.bin_peak.h5"
		"Calculate DNA methylation or chromatin accessibility epigenetic signal in spatial.")
	group_input.add_argument("--epi_feature", dest = "epi_feature", default = None,
		help = "Location of WCG/GCH/bin.merge.peak")

	group_output = workflow.add_argument_group("Output arguments")
	group_output.add_argument("--out_dir", dest = "out_dir", default = ".", 
		help = "Path to the directory where the result file shall be stored. DEFAULT: current directory. ")
	group_output.add_argument("--out_prefix", dest = "out_prefix", choices = ['WCG','GCH'], default = "WCG", 
		help = "Prefix of output files. WCG or GCH. "
		"If not specified, MultiSpace will set WCG as default.")

	group_model = workflow.add_argument_group("Model arguments")
	group_model.add_argument("--sc-scale-factor", dest = "sc_scale_factor", type = float, default = None,
		help = "The scale factor for cell-level normalization. For example, 10000. "
		"If not specified, MultiSpace will set the 75%% quantile of nCount as default. ")
	group_model.add_argument("--st-scale-factor", dest = "st_scale_factor", type = float, default = None,
		help = "The scale factor for spot-level normalization. For example, 10000. "
		"If not specified, MultiSpace will set the 75%% quantile of nCount for ST as default. ")
	group_model.add_argument("--normalize", dest = "normalize", action = "store_true", 
		help = "Whether or not to normalize the single-cell and the spatial count matrix. "
		"If set, the two matrices will be normalized by the SD for each gene. ")
	group_model.add_argument("--ntopics", dest = "ntopics_list", default = [], nargs = "+", type = int, 
		help = "Number of topics to train and test the model. MultiSpace will automatically select the optimal topic number. "
		"Multiple numbers should be separated by space. For example, --ntopics 6 7 8 9 10 . "
		"If not specified, MultiSpace will run several models with different topic numbers, and select the optimal one. ")
 

def get_mean_gene_expr(a,df):
	mean_ar = np.array(df[df.columns.intersection(a)].mean(numeric_only = True,axis = 1))
	return(mean_ar)

def get_meanEpi_methy_expr(a,df):
	mean_ar = df[df.columns.intersection(a)].mean(axis = 1).astype('float16')
	return(mean_ar)



def Mapping(sc_count_file , sc_anno_file , st_count_file ,  sc_scale_factor , st_scale_factor, out_dir , out_prefix, normalize , model_dir, gene_use , ntopics_list, spatial_location, epi_binfile, epi_feature):


	spatial_array = Deconvolve(sc_count_file, sc_anno_file, st_count_file, sc_scale_factor, st_scale_factor, out_dir, out_prefix, normalize, model_dir , gene_use , ntopics_list )

	allfiles = os.listdir(out_dir)
	peakpattern = []
	pattern =  "_topic_spot_mat_"

	for file in allfiles:
		if re.search(pattern, file):
			peakpattern.append(file)
		else:
			pass

	topic_spot_df = pd.read_csv(os.path.join(out_dir, peakpattern[0]), sep = "\t", index_col = 0)
	ntopic = peakpattern[0].split("_topic_spot_mat_")[1].split(".")[0]

	model_dir = os.path.join(out_dir, "model/")
	topic_cell_df = pd.read_csv(os.path.join(model_dir, "topic_cell_mat_%s.txt" % ntopic), sep = "\t", index_col = 0)
	cell_gene_df = pd.read_csv(sc_count_file,sep = "\t")
	spatial_count = pd.read_csv(st_count_file, sep = "\t")
	spa_loc = pd.read_csv(spatial_location, sep = "\t")


	dist_mat = distance_matrix(topic_spot_df.T, topic_cell_df.T)

	# calculate distance by topic and sort 
	dist_mat_argsort = np.argsort(dist_mat, axis = 1)
	dist_mat_argsort_top = dist_mat_argsort[:,0:40]

	dist_df_argsort_top = pd.DataFrame(dist_mat_argsort_top, index = topic_spot_df.columns)
	spot_cell_df = dist_df_argsort_top.apply(lambda x: topic_cell_df.columns[x])


	rest_topic_cell = topic_cell_df.loc[:,~topic_cell_df.columns.isin(np.unique(spot_cell_df.to_numpy().tolist()))].columns
	rest_cols_index = [topic_cell_df.columns.get_loc(col) for col in rest_topic_cell]
	cell2spot_dist_mat_argsort_top = np.argsort(dist_mat[:,rest_cols_index],axis = 0)[0:1,:]
	cell2spot_dist_mat_name_top = [topic_spot_df.columns[i] for i in cell2spot_dist_mat_argsort_top]


	rest_cell2spot_dict = defaultdict(list)
	for k,v in zip(cell2spot_dist_mat_name_top[0],rest_topic_cell):
		rest_cell2spot_dict[k].append(v)


	spot_cell_dict = spot_cell_df.T.to_dict('list')


	spot_all_cell_dict = defaultdict(list)
	for d in (spot_cell_dict, rest_cell2spot_dict): 
		for key, value in d.items():
			spot_all_cell_dict[key].append(value)

	print("Mapping Cell to Spatial...")
	spot_gene_df = []
	for key, value in spot_all_cell_dict.items():
		spot_gene_df.append(get_mean_gene_expr(value[0],cell_gene_df))
	spot_gene_df = pd.DataFrame(np.array(spot_gene_df)).T


	spot_gene_df.columns = spot_all_cell_dict.keys()
	spot_gene_df.index = cell_gene_df.index

	# spot_gene_df.to_csv(os.path.join(out_dir, "spot_gene_df.csv"))

	print("Getting Epigenetic Bin File...")
	(Epi_df,cnum) = Get_EpiSignal(epi_binfile, epi_feature, cell_gene_df)

	print("Getting Spatial Epigenetic Signal...")
	Epi_methy_df = []
	for key, value in spot_all_cell_dict.items():
	    Epi_methy_df.append(get_meanEpi_methy_expr(value[0],Epi_df))


	rm_row = np.where(~((np.asarray(Epi_methy_df).T > 0.2).sum(axis = 1) >= 10))
	Epi_methy_df = np.asarray(Epi_methy_df).T[(np.asarray(Epi_methy_df).T > 0.2).sum(axis = 1) >= 10]

	print("Saving Results...")
	sparse.save_npz(os.path.join(out_dir + out_prefix + ".signal_mat.npz"), sparse.csc_matrix(Epi_methy_df))
	pd.DataFrame(np.delete(Epi_df.index, rm_row[0])).to_csv(os.path.join(out_dir + out_prefix + ".signal_mat_rowname.txt"),index=None,header=None)


def Get_EpiSignal(epi_binfile, epi_feature, cell_gene_df):

	mat = h5py.File(epi_binfile,'r')
	mat = mat['Mcsc']
	M1 = sparse.csc_matrix((mat['data'][:],mat['indices'][:],mat['indptr'][:]), mat.attrs['shape'])
	features = pd.read_csv(epi_feature,header=None)
	Epi_df = pd.DataFrame.sparse.from_spmatrix(M1)


	Epi_df.index = features[0]
	Epi_df.columns = cell_gene_df.columns

	cnum = round(M1.shape[1]*0.05)
	Epi_df = Epi_df[(Epi_df > 0 ).sum(axis=1) >= cnum]

	return(Epi_df, cnum)








