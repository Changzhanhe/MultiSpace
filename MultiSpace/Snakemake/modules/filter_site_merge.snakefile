#filter methylation site minmum depth = 2

checkpoint read_usecells:
	input:
		config['directory']  + "QCtable.txt"
	output:
		config['directory'] + "usecells.txt"
	run:
		table = pd.read_csv(input[0],sep = "\t")
		rmlist = table[(table.GCHrownum < 500000) | (table.WCGrownum < 50000)].cell
		table[~table.isin(rmlist)['cell']]['cell'].to_csv(output[0],sep = "\t", index = False ,header = None)



def getusecells(wildcards):
	usecell = []
	with open(checkpoints.read_usecells.get(**wildcards).output[0]) as f:
		lines = f.readlines()  
		for line in lines:
			usecell.append(line.strip().replace("\n",""))
	return usecell



def site_target(wildcards):
	ls = []
	for sample in getusecells(wildcards):
		ls.append(config['dnadir'] + config['site'] + "%s/%s.WCG.sorted.merge.bed" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.WCG.binary.bed" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.GCH.sorted.merge.bed" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.GCH.chr1_chr6.binary.bed" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.GCH.chr7_chr12.binary.bed" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.GCH.chr13_chrY.binary.bed" % (sample,sample))
	ls.append(config['dnadir'] + config['site'] + "WCG.uniq.peak")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr1_chr6.uniq.peak")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr7_chr12.uniq.peak")
	ls.append(config['dnadir'] + config['site'] + "GCH.chr13_chrY.uniq.peak")
	ls.append(config['dnadir'] + config['site'] + "WCG.bin.merge.peak")
	ls.append(config['dnadir'] + config['site'] + "GCH.bin.merge.peak")
	ls.append(config['dnadir'] + config['site'] + "usecells.txt")
	return ls




rule regionbed:
	output:
		config['dnadir'] + config['site'] + "WCG.region.bed",
		config['dnadir'] + config['site'] + "GCH.region.bed"
	params:
		ref_fai=config['fasta_fai']
	shell:
		"bedtools makewindows -g {params.ref_fai} -w 500 | sort -k1,1 -k2,2n - | grep -v random - > {output[0]}  && \
		 bedtools makewindows -g {params.ref_fai} -w 1000 | sort -k1,1 -k2,2n - | grep -v random - > {output[1]}"




rule GCHmr2binary:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr1_chr6",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr7_chr12",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr13_chrY"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr1_chr6.sorted.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr7_chr12.sorted.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr13_chrY.sorted.bed",
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr1_chr6.temp"),
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr7_chr12.temp"),
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr13_chrY.temp"),
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr1_chr6.binary.bed",		
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr7_chr12.binary.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr13_chrY.binary.bed"
	params:
		binary=config['Site_binary']
	message: "TRANSFORM TO BINARY FILE: {wildcards.sample}.GCH.binary.bed"
	shell:
		"""cat {input[0]} | awk '{{OFS="\\t"}}{{if ($5 >= 2)print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n - > {output[0]}    &&\
		   cut -f 1,2,3,5 {output[0]} > {output[3]}                                                                          &&\
		   python {params.binary} --peak {output[3]} --out {output[6]}  --type GCH                                           &&\
		   cat {input[1]} | awk '{{OFS="\\t"}}{{if ($5 >= 2)print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n - > {output[1]}    &&\
		   cut -f 1,2,3,5 {output[1]} > {output[4]}                                                                          &&\
		   python {params.binary} --peak {output[4]} --out {output[7]}  --type GCH                                           &&\
		   cat {input[2]} | awk '{{OFS="\\t"}}{{if ($5 >= 2)print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n - > {output[2]}    &&\
		   cut -f 1,2,3,5 {output[2]} > {output[5]}                                                                          &&\
		   python {params.binary} --peak {output[5]} --out {output[8]}  --type GCH                                           &&\
		   rm {input[0]} &&\
		   rm {input[1]} &&\
		   rm {input[2]}"""



rule WCGmr2binary:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr1_chr6",
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr7_chr12",
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr13_chrY"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.site.sorted.bed",
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.site.temp"),
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.binary.bed"
	params:
			binary=config['Site_binary']
	shell:
		"""cat {input[0]} {input[1]} {input[2]} | awk '{{OFS="\\t"}}{{if ($5 >= 2)print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n - > {output[0]}    &&\
		   cut -f 1,2,3,5 {output[0]}  > {output[1]}                                                   &&\
		   python {params.binary} --peak {output[1]} --out {output[2]}  --type WCG                     &&\
		   rm {input[0]} &&\
		   rm {input[1]} &&\
		   rm {input[2]} """



#WCG 500bp >= 1 site
rule WCGmergesite:
	input:
		config['dnadir'] + config['site'] + "WCG.region.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.site.sorted.bed"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.sorted.merge.bed",
	shell:
		"""bedtools map -a {input[0]} -b {input[1]} -o mean | awk '{{FS="\\t";OFS="\\t"}}{{if ($4 >= 0.5){{print $1"_"$2"_"$3,1}} if ($4 != "." && $4 < 0.5){{print $1"_"$2"_"$3,0}}}}' - > {output}  """




#GCH 1000bp >= 3 site
rule GCHmergesite:
	input:
		config['dnadir'] + config['site'] + "GCH.region.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr1_chr6.sorted.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr7_chr12.sorted.bed",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.site.chr13_chrY.sorted.bed",
	output:
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.sorted.chr1_chr6.merge.bed"),		
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.sorted.chr7_chr12.merge.bed"),
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.sorted.chr13_chrY.merge.bed"),
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.sorted.merge.bed"
	shell:
		"""bedtools map -a {input[0]} -b {input[1]} -o mean,count | awk '{{FS="\\t";OFS="\\t"}}{{if ($5 >=3 && $4 >= 0.5){{print $1"_"$2"_"$3,1}} if ($5 >=3 && $4 != "." && $4 < 0.5){{print $1"_"$2"_"$3,0}}}}' - > {output[0]}  &&\
		   bedtools map -a {input[0]} -b {input[2]} -o mean,count | awk '{{FS="\\t";OFS="\\t"}}{{if ($5 >=3 && $4 >= 0.5){{print $1"_"$2"_"$3,1}} if ($5 >=3 && $4 != "." && $4 < 0.5){{print $1"_"$2"_"$3,0}}}}' - > {output[1]}  &&\
		   bedtools map -a {input[0]} -b {input[3]} -o mean,count | awk '{{FS="\\t";OFS="\\t"}}{{if ($5 >=3 && $4 >= 0.5){{print $1"_"$2"_"$3,1}} if ($5 >=3 && $4 != "." && $4 < 0.5){{print $1"_"$2"_"$3,0}}}}' - > {output[2]}  &&\
		   cat {output[0]} {output[1]} {output[2]} >> {output[3]} &&\
		   rm {input[0]} &&\
		   rm {input[1]} &&\
		   rm {input[2]}"""




#All WCG site
rule mergeWCGbinary:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.binary.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "WCG.peak"),
		temp(config['dnadir'] + config['site'] + "WCG.unsorted.peak"),
		config['dnadir'] + config['site'] + "WCG.uniq.peak"
	shell:
		"""cut -f 1 {input} >>  {output[0]}                                             &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			                            &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}                       """




rule mergeGCHbinarychr1_chr6:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr1_chr6.binary.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "GCH.chr1_chr6.peak"),
		temp(config['dnadir'] + config['site'] + "GCH.chr1_chr6.unsorted.peak"),
		config['dnadir'] + config['site'] + "GCH.chr1_chr6.uniq.peak"
	shell:
		"""cut -f 1 {input} >>  {output[0]}                                             &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			                            &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}                       """



rule mergeGCHbinarychr7_chr12:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr7_chr12.binary.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "GCH.chr7_chr12.peak"),
		temp(config['dnadir'] + config['site'] + "GCH.chr7_chr12.unsorted.peak"),
		config['dnadir'] + config['site'] + "GCH.chr7_chr12.uniq.peak"
	shell:
		"""cut -f 1 {input} >>  {output[0]}                                             &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			                            &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}                       """



rule mergeGCHbinarychr13_chrY:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr13_chrY.binary.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "GCH.chr13_chrY.peak"),
		temp(config['dnadir'] + config['site'] + "GCH.chr13_chrY.unsorted.peak"),
		config['dnadir'] + config['site'] + "GCH.chr13_chrY.uniq.peak"
	shell:
		"""cut -f 1 {input} >>  {output[0]}                                             &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			                            &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}                       """





rule WCGmergebin:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.sorted.merge.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "WCG.bin.merge.temp"),
		temp(config['dnadir'] + config['site'] + "WCG.bin.merge.unsort.peak"),
		config['dnadir'] + config['site'] + "WCG.bin.merge.peak"
	shell:
		"""cut -f 1 {input} >> {output[0]}                           &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			         &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}    """




rule GCHmergebin:
	input:
		lambda wildcards: expand(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.sorted.merge.bed", sample = getusecells(wildcards))
	output:
		temp(config['dnadir'] + config['site'] + "GCH.bin.merge.temp"),
		temp(config['dnadir'] + config['site'] + "GCH.bin.merge.unsort.peak"),
		config['dnadir'] + config['site'] + "GCH.bin.merge.peak"
	shell:
		"""cut -f 1 {input} >> {output[0]}                           &&\
		awk '!a[$0]++' {output[0]}	> {output[1]}			         &&\
		cat {output[1]} | sort -t "_" -k1,1 -k2,2n  > {output[2]}    """




# rule bed2bw:
# 	input:
# 		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.binary.bed",
# 		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.binary.bed"
# 	output:
# 		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.bedGraph"),
# 		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.bedGraph"),
# 		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.bw",
# 		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.bw"
# 	params:
# 		bedGraphToBigWig=config['bedGraphToBigWig'],
# 		ref_fai=config['fasta.fai']
# 	message: "CHANGE BEDGRAPH TO BIGWIG: {wildcards.sample} bw file"
# 	shell:
# 		"cut -f 1,2,3,5 {input[0]} > {output[0]}                             &&\
# 		{params.bedGraphToBigWig} {output[0]} {params.ref_fai} {output[2]}   &&\
# 		 cut -f 1,2,3,5 {input[1]} > {output[1]}                             &&\
# 		{params.bedGraphToBigWig} {output[1]} {params.ref_fai} {output[3]}    "



