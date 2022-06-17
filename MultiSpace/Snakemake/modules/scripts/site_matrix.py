#!/usr/bin/env python

import os,sys
import pandas as pd
import h5py
import numpy as np
import scipy.sparse as sp_sparse

from scipy import sparse,io
from scipy.sparse import csr_matrix
from pandas import Series,DataFrame
from optparse import OptionParser



def prepare_optparser():
	usage = "USAGE: %prog -c 04.WCG.GCH/usecells.txt -p 04.WCG.GCH/ -t [GCH/WCG]"
	optparser = OptionParser(usage=usage)
	optparser.add_option("-c", "--usecell", help="path to usecell.txt")
	optparser.add_option("-p", "--sitepath", help="path to uniq.peak")
	optparser.add_option("-t", "--type", help="WCG or GCH")
	optparser.add_option("-r", "--chr", help="chrX1_chrX2")
	return optparser


def main():
	parser = prepare_optparser()
	(options,args) = parser.parse_args()

	try:
		barcode_list = pd.read_csv(options.usecell,header = None)
		sitepath = options.sitepath
		sitetype = options.type
		chrname = options.chr

	except IndexError:
		prepare_optparser().print_help()
		sys.exit(1)

	if sitetype == "WCG" :
		WCGsite(barcode_list, sitepath, sitetype)
	elif sitetype == "GCH":
		GCHsite(barcode_list, sitepath, sitetype, chrname)



def GCHsite(barcode_list, sitepath, sitetype, chrname):

	site_peak = pd.read_table(os.path.join(sitepath + sitetype + "." + chrname +".uniq.peak"), header = None, names = ["name"], index_col = False)

	size = 0
	indptr  = [0]
	indices    = []
	data = []


	for index, row in barcode_list.iterrows():
		tmp = pd.read_table(sitepath + "%s/%s.%s.%s.binary.bed" % (row[0], row[0], sitetype, chrname), names = ["name","met"], index_col = False)
		print( "%s/%s.%s.%s.binary.bed" % (row[0],row[0], sitetype, chrname))
		tmp = pd.concat([site_peak.set_index('name'),tmp.set_index('name')],axis = 1,join = "outer").reset_index()
		size += tmp.loc[~tmp['met'].isnull()].iloc[:,0].size
		indptr = np.append(indptr,size)
		indices = np.append(indices,tmp.loc[~tmp['met'].isnull()].index)
		data = np.append(data,np.array(tmp[~tmp['met'].isnull()]['met']))


	write_h5(filename = os.path.join(sitepath + sitetype + "." + chrname +".site_peak.h5"), data = data, indptr = indptr, indices = indices, peakmat = site_peak)


def WCGsite(barcode_list, sitepath, sitetype):

	site_peak = pd.read_table(os.path.join(sitepath + sitetype + ".uniq.peak"), header = None, names = ["name"], index_col = False)

	size = 0
	indptr  = [0]
	indices    = []
	data = []


	for index, row in barcode_list.iterrows():
		tmp = pd.read_table(sitepath + "%s/%s.%s.binary.bed" % (row[0], row[0], sitetype), names = ["name","met"], index_col = False)
		print("%s/%s.%s.binary.bed" % (row[0],row[0], sitetype))
		tmp = pd.concat([site_peak.set_index('name'),tmp.set_index('name')],axis = 1,join = "outer").reset_index()
		size += tmp.loc[~tmp['met'].isnull()].iloc[:,0].size
		indptr = np.append(indptr,size)
		indices = np.append(indices,tmp.loc[~tmp['met'].isnull()].index)
		data = np.append(data,np.array(tmp[~tmp['met'].isnull()]['met']))


	write_h5(filename = os.path.join(sitepath + sitetype +".site_peak.h5"), data = data, indptr = indptr, indices = indices, peakmat = site_peak)



def write_h5(filename, data, indptr, indices, peakmat):

	f = h5py.File(filename,'w')
	mat = f.create_group('Mcsc')
	mat.create_dataset('data',data=data)
	mat.create_dataset('indptr',data=indptr)
	mat.create_dataset('indices',data=indices)
	mat.attrs['shape'] = (len(peakmat),len(indptr)-1)
	f.close()


if __name__=='__main__':
	main()





