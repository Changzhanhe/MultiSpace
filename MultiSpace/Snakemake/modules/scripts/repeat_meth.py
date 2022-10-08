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
	optparser.add_option("-e", "--repeat", choices = ['methylation', 'expression'], help = "methylation or expression")
	return optparser


def main():
	parser = prepare_optparser()
	(options,args) = parser.parse_args()

	try:
		barcode_list = pd.read_csv(options.usecell,header = None)
		path = options.path
		methtype = options.methtype
		repeat = options.repeat

	except IndexError:
		prepare_optparser().print_help()
		sys.exit(1)

	if repeat == "methylation" :
		meth(barcode_list, path, methtype)
	elif repeat == "expression":
		expr(barcode_list, path, methtype)



def meth(barcode_list, path, methtype):

	for index, row in barcode_list.iterrows():
		file = pd.read_table(os.path.join(path + "%s/%s.repeat.%s" % (row[0], row[0], methtype)),header = None)
		file.columns = ["name","met"]
		listdf = pd.DataFrame(pd.concat([listdf.set_index('name'),file.set_index('name')],axis = 1,join = "outer").reset_index())

	listdf.index = listdf["name"]
	listdf.drop('name', inplace=True, axis=1)

	listdf.to_csv(os.path.join(path, "%smeth.csv" % (methtype)))

def expr(barcode_list, path, methtype):

	listdf = pd.DataFrame(columns=["ID"])

	for index, row in barcode_list.iterrows():
		file = pd.read_table(os.path.join(path + "%s.repeats.%s.txt" % (row[0], methtype)),index_col = False).iloc[:,[0,8]]
		file.columns = ['ID',row[0]]
		listdf = pd.DataFrame(pd.concat([listdf.set_index('ID'),file.set_index('ID')],axis = 1,join = "outer").reset_index())

	listdf.index = listdf['ID']
	listdf.drop('ID', axis = 1 ,inplace = True)
	listdf.index = listdf.index.str.split('|',1).str[0]
	listdf = listdf.apply(lambda x: (x/x.sum())*1e6, axis = 0)

	listdf.to_csv(os.path.join(path, "%sexpr_norm.csv" % (methtype)))


if __name__=='__main__':
	main()


