#calling methylation site minmum depth = 1
def methy_target(wildcards):
	ls = []
	for sample in config['runs']:
		# ls.append(config['dnadir'] + config['site'] + "%s/%s.CG.mr" % (sample,sample))
		ls.append(config['dnadir'] + config['site'] + "%s/%s.QC.txt" % (sample,sample))
	ls.append(config['dnadir'] + config['site'] + "QCtable.txt")
	ls.append(config['dnadir'] + config['site'] + "usecells.txt")	
	return ls


def get_rmdup_sort_bam(wildcards):
	bam = []
	bam.append(config['dnadir'] + config['bam'] + "%s/%s.sort.rmdup.bam" % (wildcards.sample,wildcards.sample))
	return bam



rule methylation_call_chr1_chr6:
	input: 
		get_rmdup_sort_bam,
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr1_chr6"
	params:
		ref_fa=config['lambda_fasta'],
		scripts = config['methratio']
	message: "CALL METH: {wildcards.sample} methylation calling: chr1 - chr6"
	shell:
		"python2 {params.scripts} -o {output} -z -n -q --chr=chr1,chr2,chr3,chr4,chr5,chr6 -d {params.ref_fa} {input}"



rule methylation_call_chr7_chr12:
	input: 
		get_rmdup_sort_bam,
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr7_chr12"
	params:
		ref_fa=config['lambda_fasta'],
		scripts = config['methratio']
	message: "CALL METH: {wildcards.sample} methylation calling: chr7 - chr12"
	shell:
		"python2 {params.scripts} -o {output} -z -n -q --chr=chr7,chr8,chr9,chr10,chr11,chr12 -d {params.ref_fa} {input}"



rule methylation_call_chr13_chrY:
	input: 
		get_rmdup_sort_bam,
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr13_chrY"
	params:
		ref_fa=config['lambda_fasta'],
		scripts = config['methratio'] 
	message: "CALL METH: {wildcards.sample} methylation calling: chr13 - chrY"
	shell:
		"python2 {params.scripts} -o {output} -z -n -q --chr=chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrM,chrX,chrY -d {params.ref_fa} {input}"



rule conversion_rate:
	input:
		get_rmdup_sort_bam,
	output:
		temp(config['dnadir'] + config['site'] + "{sample}/{sample}.CG.lambda"),
		config['dnadir'] + config['site'] + "{sample}/{sample}.QC.txt"
	params:
		ref_fa=config['lambda_fasta'],
		scripts = config['methratio']
	message: "CALL METH: {wildcards.sample} calculate conversion rate"
	shell:
		"""python2 {params.scripts} -o {output[0]} -z -n -q --chr=lambda -d {params.ref_fa} {input}   &&\
		cat {output[0]} | awk '{{if($3=="+" && $4 ~/[AGCT][AT]CG[AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT]CG[AT][AGCT]/){{print $0}}}}'    |\
		awk 'BEGIN{{sum1=0;sum2=0}}{{sum1+=$6;sum2+=$7}}END{{print "'{wildcards.sample}'\\tWCG\\t"sum2/sum1}}' >> {output[1]}                         &&\
		cat {output[0]} | awk '{{if($3=="+" && $4 ~/[AGCT]GC[ACT][AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT][AGT]GC[AGCT]/){{print $0}}}}'      |\
		awk 'BEGIN{{sum1=0;sum2=0}}{{sum1+=$6;sum2+=$7}}END{{print "'{wildcards.sample}'\\tGCH\\t"sum2/sum1}}' >> {output[1]}                         """




rule get_methysite_chr1_chr6:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr1_chr6"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr1_chr6",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr1_chr6",
	message: "GET METHYLATION SITE: {wildcards.sample} WCG and GCH: chr1 - chr6"
	shell:
		"""cat {input} | grep -v lambda | grep -v random | awk '{{if ($3=="+" && $4 ~/[AGCT][AT]CG[AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT]CG[AT][AGCT]/){{print $0}}}}' |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[0]}       &&\
		cat {input} | grep -v lambda | grep -v random | awk '{{if($3=="+" && $4 ~/[AGCT]GC[ACT][AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT][AGT]GC[AGCT]/){{print $0}}}}'   |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[1]}   &&\
		rm {input}"""



rule get_methysite_chr7_chr12:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr7_chr12"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr7_chr12",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr7_chr12",
	message: "GET METHYLATION SITE: {wildcards.sample} WCG and GCH: chr7 - chr12"
	shell:
		"""cat {input} | grep -v lambda | grep -v random | awk '{{if ($3=="+" && $4 ~/[AGCT][AT]CG[AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT]CG[AT][AGCT]/){{print $0}}}}' |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[0]}       &&\
		cat {input} | grep -v lambda | grep -v random | awk '{{if($3=="+" && $4 ~/[AGCT]GC[ACT][AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT][AGT]GC[AGCT]/){{print $0}}}}'   |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[1]}   &&\
		rm {input}"""



rule get_methysite_chr13_chrY:
	input:
		config['dnadir'] + config['site'] + "{sample}/{sample}.CG.chr13_chrY"
	output:
		config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.chr13_chrY",
		config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.chr13_chrY",
	message: "GET METHYLATION SITE: {wildcards.sample} WCG and GCH: chr13 - chrY"
	shell:
		"""cat {input} | grep -v lambda | grep -v random | awk '{{if ($3=="+" && $4 ~/[AGCT][AT]CG[AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT]CG[AT][AGCT]/){{print $0}}}}' |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[0]}       &&\
		cat {input} | grep -v lambda | grep -v random | awk '{{if($3=="+" && $4 ~/[AGCT]GC[ACT][AGCT]/){{print $0}}else if($3=="-" && $4 ~/[AGCT][AGT]GC[AGCT]/){{print $0}}}}'   |\
		cut -f 1,2,3,5,6,7 - | grep -v NA > {output[1]}   &&\
		rm {input}"""



rule GCH_site_num:
	input:
		expand(config['dnadir'] + config['site'] + "{sample}/{sample}.GCH.{ext}", sample = config['runs'],ext = ['chr1_chr6','chr7_chr12','chr13_chrY']),
	output:
		config['dnadir'] + config['site'] + "GCHsite.num"
	shell:
		"wc -l {input} >> {output}  "



rule WCG_site_num:
	input:
		expand(config['dnadir'] + config['site'] + "{sample}/{sample}.WCG.{ext}", sample = config['runs'],ext = ['chr1_chr6','chr7_chr12','chr13_chrY']),
	output:
		config['dnadir'] + config['site'] + "WCGsite.num"
	shell:
		"wc -l {input} >> {output}  "



# rule mergeQCandfiltercell:
# 	input:
# 		expand(config['dnadir'] + config['site'] + "{sample}/{sample}.QC.txt",  sample = config['runs'])
# 	output:
# 		config['dnadir'] + config['site'] + "sample.Conv.txt",
# 		config['dnadir'] + config['site'] + "QCtable.txt",
# 		config['dnadir'] + config['site'] + "usecells.txt"	
# 	params:
# 		scripts = config['QCtable'],
# 	shell:
# 		"cat {input} >> {output[0]}  && python {params.scripts}"



rule mergeQC:
	input:
		expand(config['dnadir'] + config['site'] + "{sample}/{sample}.QC.txt",  sample = config['runs'])
	output:
		config['dnadir'] + config['site'] + "sample.Conv.txt"
	shell:
		"cat {input} >> {output}"



rule mergeQCandfiltercell:
	input:
		config['dnadir'] + config['site'] + "WCGsite.num",
		config['dnadir'] + config['site'] + "GCHsite.num",
		config['dnadir'] + config['site'] + "sample.Conv.txt"
	output:
		config['directory']  + "QCtable.txt",
		# config['dnadir'] + config['site'] + "usecells.txt"	
	params:
		scripts = config['QCtable'],
	shell:
		"python {params.scripts} -W {input[0]} -G {input[1]} -Q {output}  &&\
		rm {input[0]} &&\
		rm {input[1]} &&\
		rm {input[2]}"


rule rmbam:
	input:
		config['dnadir'] + config['bam']
	shell:
		"rm -rf {input}"

