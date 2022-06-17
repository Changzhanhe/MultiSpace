# -*- coding: utf-8 -*-
# @Author: Zhanhe Chang

import os,sys
import re
import pandas as pd
import shutil
from optparse import OptionParser




def main():
	usage = "USAGE: %prog -W [path to WCGsite.num] -G [path to GCHsite.num]"
	optparser = OptionParser(usage=usage)
	optparser.add_option("-W", "--WCG", help="path to WCGsite.num")
	optparser.add_option("-G", "--GCH", help="path to GCHsite.num")
	optparser.add_option("-Q", "--QC", help="path to output QCtable")
	# optparser.add_option("-C", "--CELL", help="path to usecells")
	(options, args) = optparser.parse_args(sys.argv)



	bamworkpath = []
	siteworkpath = []
	with open("./config.yaml") as f:
		for line in f:
			if line.strip().split(":")[0] == "coolseqdir":
				bamworkpath.append(line.strip().split(":")[1].replace("\"","").replace(" ","") + "03.Bam/")  
				siteworkpath.append(line.strip().split(":")[1].replace("\"","").replace(" ","") + "04.WCG.GCH/")



	listOfFiles = []
	lst = []
	ra = []
	mean = []

	for dirpath, dirnames, filenames in os.walk(bamworkpath[0]):
		listOfFiles = [os.path.join(dirpath,file) for file in filenames]
		for log in listOfFiles:
			if log.endswith(".log"):
				lst.append(re.split('\/',log)[-2])
				f = open(log,'r')
				lines = f.readlines()  
				for line in lines:
					if "aligned reads" in line:
						ra.append(float(re.split(':| |\(',line)[5].replace("%","").replace(")","").replace(",","")))


	step = 2
	ra_ = [ra[i:i+step] for i in range(0,len(ra), step)]
	for num in ra_:
		mean.append(round(num[0] + num[1])/2)


	d = {'cell':lst, 'mean_map_ratio':mean}
	df = pd.DataFrame(d, columns=["cell","mean_map_ratio"])


	QC = pd.read_csv(siteworkpath[0] + "sample.Conv.txt", sep = "\t", header = None, names = ['cell','type','con'])


	group_obj = QC.groupby('type')
	for key, item in group_obj:
	    df = pd.merge(pd.DataFrame(item), df, on = "cell")
	df = df[['cell','con_x','con_y','mean_map_ratio']]


	df.columns = ['cell','WCG_conv','GCH_conv','mean_map_ratio']


	list = []
	with open(options.WCG) as site:
		lines = site.readlines()
		lines = lines[:-1]
		for line in lines:
			dict = {}
			dict['cell'] = line.strip().split("/")[-2]
			dict['WCGrownum'] = int(line.strip().split("/")[0].replace(" ",""))
			list.append(dict)


	table = pd.merge(pd.DataFrame(list).groupby('cell').sum(),df, on = "cell")

	list = []
	with open(options.GCH) as site:
		lines = site.readlines()
		lines = lines[:-1]
		for line in lines:
			dict = {}
			dict['cell'] = line.strip().split("/")[-2]
			dict['GCHrownum'] = int(line.strip().split("/")[0].replace(" ",""))
			list.append(dict)


	table = pd.merge(pd.DataFrame(list).groupby('cell').sum(),table, on = "cell")

	table.to_csv(os.path.join(options.QC), sep = "\t", index = False )

	# rmlist = table[(table.GCHrownum < 500000) | (table.WCGrownum < 50000)].cell
	# table[~table.isin(rmlist)['cell']]['cell'].to_csv(os.path.join(siteworkpath[0] + "usecells.txt"),sep = "\t", index = False ,header = None)


if __name__=='__main__':
	main()




