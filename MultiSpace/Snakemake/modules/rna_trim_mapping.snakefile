checkpoint read_rnausecells:
	input:
		config['directory']  + "04.WCG.GCH/QCtable.txt"
	output:
		config['directory'] + "04.WCG.GCH/usecells.txt"
	run:
		table = pd.read_csv(input[0],sep = "\t")
		rmlist = table[(table.GCHrownum < 500000) | (table.WCGrownum < 50000)].cell
		table[~table.isin(rmlist)['cell']]['cell'].to_csv(output[0],sep = "\t", index = False ,header = None)

def get_dnausecells(wildcards):
	usecells = []
	with open(checkpoints.read_rnausecells.get(**wildcards).output[0]) as f:
		lines = f.readlines()  
		for line in lines:
			usecells.append(line.strip().replace("\n",""))
	return usecells


def rna_target(wildcards):
	ls = []
	for sample in get_dnausecells(wildcards):
		ls.append(config['rnadir'] + config['fastq'] + "%s/%s_val_1.fq.gz" % (sample,sample)),
		ls.append(config['rnadir'] + config['fastq'] + "%s/%s_val_2.fq.gz" % (sample,sample)),
		ls.append(config['rnadir'] + config['bam'] + "%sAligned.sortedByCoord.out.bam" % (sample)),
		ls.append(config['dnadir'] + config['site'] + "LINEexpr_norm.csv"),
		ls.append(config['dnadir'] + config['site'] + "LTRexpr_norm.csv"),
		ls.append(config['dnadir'] + config['site'] + "SINEexpr_norm.csv"),
		ls.append(config['rnadir'] + config['expr'] + "rawcount.txt")
	return ls



def get_fastq(wildcards):
	return glob.glob(config['rnadir'] + config['raw'] + wildcards.sample + "_[1-2].fastq.gz")


rule rna_trim_adapter:
	input: get_fastq
	output:
		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_1.fq.gz",
		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_2.fq.gz",
		temp(config['rnadir'] + config['fastq'] + "{sample}/{sample}_1.fastq.gz_trimming_report.txt"),
		temp(config['rnadir'] + config['fastq'] + "{sample}/{sample}_2.fastq.gz_trimming_report.txt")
	params:
		basename=lambda wildcards: "%s" % (wildcards.sample),
		quality=20,
		stringency=3,
		length=75,
		output_dir=lambda wildcards: config['rnadir'] + config['fastq'] + "%s/" % (wildcards.sample)
	message: "RNA TRIM: Trim adaptors for paired-end using trim_galore - {wildcards.sample}"
	shell:
		"trim_galore --trim-n -quality {params.quality} --length {params.length} --stringency {params.stringency} --basename {params.basename} -j 7 --paired --gzip -o {params.output_dir} {input} "



rule star_mapping:
	input: 
		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_1.fq.gz",
		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_2.fq.gz",
	output:
		config['rnadir'] + config['bam'] + "{sample}Aligned.sortedByCoord.out.bam",
		temp(config['rnadir'] + config['bam'] + "{sample}Log.progress.out"),
		temp(config['rnadir'] + config['bam'] + "{sample}Log.out"),
		temp(config['rnadir'] + config['bam'] + "{sample}Log.final.out"),
		temp(config['rnadir'] + config['bam'] + "{sample}SJ.out.tab")
	params:
		ref=config['star_index'],
		thread=5,
		prefix=config['rnadir'] + config['bam'] + "{sample}",
		Nmax=999,
		Lmax=0.04,
		Mmax=500,
	message: "RNA ALIGN: Align {wildcards.sample} to the genome by STAR"
	shell:
		"""STAR --runThreadN {params.thread} --runMode alignReads --genomeDir {params.ref} --readFilesCommand zcat --readFilesIn {input[0]} {input[1]} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
		--outFilterMismatchNmax {params.Nmax} --outFilterMismatchNoverLmax {params.Lmax} --outFilterMultimapNmax {params.Mmax}"""



rule feature_count:
	input: 	
		lambda wildcards: expand(config['rnadir'] + config['bam'] + "{sample}Aligned.sortedByCoord.out.bam", sample = get_dnausecells(wildcards))
	output: 
		temp(config['rnadir'] + config['expr'] + "count.txt"),
		temp(config['rnadir'] + config['expr'] + "count.txt.summary"),
		config['rnadir'] + config['expr'] + "rawcount.txt",
		# temp(config['rnadir'] + config['bam'] +  "{sample}.log")
	params:
		featuretype="exon",
		ref_gtf=config['rna_annotation'],
		threads=5,
		attributetype="gene_id"
	message: "QUANTIFICATION: RNA quantification using FeatureCount"
	shell:
		"""featureCounts -T {params.threads} -p -t {params.featuretype} -g {params.attributetype} -a {params.ref_gtf} -o {output[0]} {input}  &&\
		cut -f 1,7- {output[0]} > {output[2]} &&\
		sed -i '1d' {output[2]}"""


rule bam2sam:
	input:
		config['rnadir'] + config['bam'] + "{sample}Aligned.sortedByCoord.out.bam",
	output:
		temp(config['rnadir'] + config['expr'] + "{sample}header.sam"),
		temp(config['rnadir'] + config['expr'] + "{sample}_repeatAligned.sortedByCoord.out.sam"),
		temp(directory(config['rnadir'] + config['expr'] + "{sample}_tags/")),
		config['rnadir'] + config['expr'] + "{sample}.repeats.LINE.txt",
		config['rnadir'] + config['expr'] + "{sample}.repeats.LTR.txt",
		config['rnadir'] + config['expr'] + "{sample}.repeats.SINE.txt",
		# temp(config['rnadir'] + config['expr'] +  "{sample}.sam.log")
	params: 
		species = config['species']
	shell:
		"""samtools view -H {input} | grep -v -w "@PG"  | grep -v -w "@CO" > {output[0]} &&\
		samtools reheader {output[0]} {input} | samtools view -h > {output[1]}  &&\
		makeTagDirectory {output[2]} {output[1]} -format sam -keepOne &&\
		analyzeRepeats.pl repeats {params.species} -d {output[2]} -L2 LINE -noadj > {output[3]} &&\
		analyzeRepeats.pl repeats {params.species} -d {output[2]} -L2 LTR -noadj > {output[4]}  &&\
		analyzeRepeats.pl repeats {params.species} -d {output[2]} -L2 SINE -noadj > {output[5]}  """



rule LINErepeatmatrix:
	input:
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.LTR.txt", sample = get_dnausecells(wildcards))
	output:
		config['rnadir'] + config['expr'] + "LTRexpr_norm.csv",
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype LTR --repeat expression"""



rule LTRrepeatmatrix:
	input:
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.LINE.txt", sample = get_dnausecells(wildcards))
	output:
		config['dnadir'] + config['site'] + "LINEexpr_norm.csv",
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype LINE --repeat expression"""


rule SINErepeatmatrix:
	input:
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.SINE.txt", sample = get_dnausecells(wildcards))
	output:
		config['dnadir'] + config['site'] + "SINEexpr_norm.csv"
	params:
		repeatmeth=config['Repeat_meth'],
		cell=config['dnadir'] + config['site'] + "usecells.txt",
		path=config['dnadir'] + config['site']
	message: "GENERATE WCG BIN BY CELL MATRIX "
	shell:
		"""python {params.repeatmeth} --usecell {params.cell} --path {params.path} --methtype SINE --repeat expression"""



rule TEGene_normalize:
	input:
		


rule rm:
	input:
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.LINE.txt", sample = get_dnausecells(wildcards)),
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.LTR.txt", sample = get_dnausecells(wildcards)),
		lambda wildcards: expand(config['rnadir'] + config['expr'] + "{sample}.repeats.SINE.txt", sample = get_dnausecells(wildcards))
	shell:
		"""rm {input[0]} && rm{input[1]}"""





# installing & updating homer software and genomes
# perl ~/biotools/homer/configureHomer.pl -install homer
# perl ~/biotools/homer/configureHomer.pl -install hg38

# rule add_umi:
# 	input: 
# 		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_1.fq.gz",
# 		config['rnadir'] + config['fastq'] + "{sample}/{sample}_val_2.fq.gz"
# 	output:
# 		config['rnadir'] + config['fastq'] + "{sample}/{sample}_R1_extracted.fq.gz",
# 		config['rnadir'] + config['fastq'] + "{sample}/{sample}_R2_extracted.fq.gz",
# 		config['rnadir'] + config['fastq'] + "{sample}/{sample}.extract.log"
# 	message: "ADD UMI: Add umi using umi_tools - {wildcards.sample}"
# 	shell:
# 		"umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN --stdin {input[0]} --read2-in {input[1]} --stdout {output[0]} --read2-out {output[1]} -L {output[2]}" 

		
#--filter-cell-barcode  
#--whitelist=$barcode                     

# if [ $type == 'single' ];then
# 	trim_galore --trim-n --length 35 -j 4 -q 20 -a AAGCAGTGGTATCAACGCAGAGTACATGGG \
#   -a2 TTTTTTTTTT --fastqc $fastq_dir/$samp.fq --gzip -o $trim_dir





