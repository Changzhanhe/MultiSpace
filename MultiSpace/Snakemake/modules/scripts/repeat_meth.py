#!/usr/bin/env python

import os,sys
import pandas as pd
import numpy as np
import functools

from pandas import Series,DataFrame
from optparse import OptionParser
from functools import reduce



def prepare_optparser():
	usage = "USAGE: %prog -c 04.WCG.GCH/usecells.txt -p 04.WCG.GCH/ -t [GCH/WCG]"
	optparser = OptionParser(usage=usage)
	optparser.add_option("-c", "--usecell", help="path to usecell.txt")
	optparser.add_option("-p", "--path", help="path to each cell folder located in")
	optparser.add_option("-t", "--methtype", help="LINE or LTR")
	return optparser



def main():
	parser = prepare_optparser()
	(options,args) = parser.parse_args()

	try:
		barcode_list = pd.read_csv(options.usecell,header = None)
		path = options.path
		methtype = options.methtype

	except IndexError:
		prepare_optparser().print_help()
		sys.exit(1)

	if methtype == "LINE" :
		LINEmeth(barcode_list, path, methtype)
	elif methtype == "LTR":
		LTRmeth(barcode_list, path, methtype)



def LINEmeth(barcode_list, path, methtype):

	list = []
	for index, row in barcode_list.iterrows():
	    list.append(os.path.join(sitepath + "%s/%s.repeat.LINE" % (row[0], row[0], methtype)))

	list_stacked = pd.DataFrame(columns=["name","met"])

	for filename in list:
	    file = pd.read_table(filename,index_col = False, names = ["name",filename])
	    list_stacked = pd.DataFrame(pd.concat([list_stacked.set_index('name'),file.set_index('name')],axis = 1,join = "outer").reset_index())

	list_stacked.index = list_stacked["name"]
	list_stacked.drop(list_stacked.iloc[:, 0:2], inplace=True, axis=1)

	list_stacked.to_csv(os.path.join(path, "LINEmeth.csv"), sep = "\t", index = False )



def LTRmeth(barcode_list, path, methtype)

	list = []
	for index, row in barcode_list.iterrows():
	    list.append(os.path.join(sitepath + "%s/%s.repeat.LTR" % (row[0], row[0], methtype)))

	list_stacked = pd.DataFrame(columns=["name","met"])

	for filename in list:
	    file = pd.read_table(filename,index_col = False, names = ["name",filename])
	    list_stacked = pd.DataFrame(pd.concat([list_stacked.set_index('name'),file.set_index('name')],axis = 1,join = "outer").reset_index())

	list_stacked.index = list_stacked["name"]
	list_stacked.drop(list_stacked.iloc[:, 0:2], inplace=True, axis=1)

	list_stacked.to_csv(os.path.join(path, "LTRmeth.csv"), sep = "\t", index = False )



if __name__=='__main__':
	main()


