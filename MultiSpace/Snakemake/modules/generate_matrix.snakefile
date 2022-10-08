
def matrix_target(wildcards):
	ls = []
	ls.append(config['dnadir'] + config['site'] + "WCG.bin_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "GCH.bin_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "WCG.site_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr1_chr6.site_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr7_chr12.site_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr13_chrY.site_peak.h5")
	ls.append(config['dnadir'] + config['site'] + "LINEmeth.csv")
	ls.append(config['dnadir'] + config['site'] + "LTRmeth.csv")
	ls.append(config['dnadir'] + config['site'] + "SINEmeth.csv")
	return ls



rule WCGbinbycell:
	input:
		config['dnadir'] + config['site'] + "WCG.bin.merge.peak"
	output:
		config['dnadir'] + config['site'] + "WCG.bin_peak.h5",
	params:
		binmat=config['Bin_matrix'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		binpath=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.binmat} --usecell {params.cell} --binpath {params.binpath} --type WCG """




rule GCHbinbycell:
	input:
		config['dnadir'] + config['site'] + "GCH.bin.merge.peak"
	output:
		config['dnadir'] + config['site'] + "GCH.bin_peak.h5",
	params:
		binmat=config['Bin_matrix'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		binpath=config['dnadir'] + config['site']
	message: "GENERATE GCH BIN BY CELL MATRIX "
	shell:
		"""python {params.binmat} --usecell {params.cell} --binpath {params.binpath} --type GCH """



rule WCGsitebycell:
	input:
		config['dnadir'] + config['site'] + "WCG.uniq.peak"
	output:
		config['dnadir'] + config['site'] + "WCG.site_peak.h5",
	params:
		sitemat=config['Site_matrix'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		sitepath=config['dnadir'] + config['site']
	message: "GENERATE WCG SITE BY CELL MATRIX "
	shell:
		"""python {params.sitemat} --usecell {params.cell} --path {params.sitepath} --type WCG --chr all"""

	

rule GCHsitebycell:
	input:
		config['dnadir'] + config['site'] + "GCH.chr1_chr6.uniq.peak",
		config['dnadir'] + config['site'] + "GCH.chr7_chr12.uniq.peak",
		config['dnadir'] + config['site'] + "GCH.chr13_chrY.uniq.peak",
	output:
		config['dnadir'] + config['site'] + "GCH.chr1_chr6.site_peak.h5",
		config['dnadir'] + config['site'] + "GCH.chr7_chr12.site_peak.h5",
		config['dnadir'] + config['site'] + "GCH.chr13_chrY.site_peak.h5",
	params:
		sitemat=config['Site_matrix'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		sitepath=config['dnadir'] + config['site']
	message: "GENERATE GCH SITE BY CELL MATRIX "
	shell:
		"""python {params.sitemat} --usecell {params.cell} --sitepath {params.sitepath} --type GCH --chr chr1_chr6 &&\
		   python {params.sitemat} --usecell {params.cell} --sitepath {params.sitepath} --type GCH --chr chr7_chr12 &&\
		   python {params.sitemat} --usecell {params.cell} --sitepath {params.sitepath} --type GCH --chr chr13_chrY """



rule LINEmethmatrix:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.repeat.LINE",
	output:
		config['dnadir'] + config['site'] + "LINEmeth.csv",
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype LINE --repeat methylation"""



rule LTRmethmatrix:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.repeat.LTR",
	output:
		config['dnadir'] + config['site'] + "LTRmeth.csv",
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype LTR --repeat methylation"""


rule SINEmethmatrix:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.repeat.SINE",
	output:
		config['dnadir'] + config['site'] + "SINEmeth.csv",
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype SINE --repeat methylation"""



rule rmsite:
	input:
		config['dnadir'] + config['site'] + "{sample}"
	shell:
		"rm {input}"




