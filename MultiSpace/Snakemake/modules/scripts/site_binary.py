import argparse
from argparse import ArgumentParser
import sys
import math


if __name__ == "__main__" :

	parser = argparse.ArgumentParser(description='Binarize peak methylation ratio matrix for each cell')
	parser.add_argument("--peak",dest = "peak", help='Peak count matrix', required = True)
	parser.add_argument("--out",dest = "out", help = 'Binarized peak',required = True)
	parser.add_argument("--type",dest = "type",help='Which type of data', required =True)
	options = parser.parse_args()



def binary_peak_count(options):

	if options.type == "WCG":
		fhd = open(options.peak, 'rt')
		with open(options.out , 'w') as f:
			for line in fhd.readlines():
				a,b,c,d = [x.strip() for x in line.split('\t')]
				if (float(d) % 1) > 0.5:
					d = math.ceil(float(d))
				if (float(d) % 1) < 0.5:
					d = round(float(d))
				else:
					d = 0.5
				f.write(str(a) + '_' + str(b) +'_' + str(c) +'\t' + str(d) +'\t' + '\n')
			f.close()
	if options.type == "GCH":
		fhd = open(options.peak,'rt')
		with open(options.out , 'w') as f:
			for line in fhd.readlines():
				a,b,c,d = [x.strip() for x in line.split('\t')]
				if (float(d) % 1) >= 0.5:
					d = math.ceil(float(d))
				else:
					d = round(float(d))
				f.write(str(a) + '_' + str(b) +'_' + str(c) +'\t' + str(d) +'\t' + '\n')
			f.close()

binary_peak_count(options)
