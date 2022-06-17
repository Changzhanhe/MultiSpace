
def get_dnausecells(wildcards):
	usecells = []
	with open(checkpoints.read_usecells.get(**wildcards).output[0]) as f:
		lines = f.readlines()  
		for line in lines:
			usecells.append(line.strip().replace("\n",""))
	return usecells


def rna_target(wildcards):
	ls = []
	for sample in get_dnausecells(wildcards):
		ls.append(config['rnadir'] + config['fastq'] + "%s/%s_val_1.fq.gz" % (sample,sample)),
		ls.append(config['rnadir'] + config['fastq'] + "%s/%s_val_2.fq.gz" % (sample,sample)),
		ls.append(config['rnadir'] + config['bam'] + "%sAligned.sortedByCoord.out.bam" % (sample))
	ls.append(config['rnadir'] + config['bam'] + "rawcount.txt")
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
	log: config['rnadir'] + config['fastq'] +  "log/{sample}.log"
	shell:
		"trim_galore --trim-n -quality {params.quality} --length {params.length} --stringency {params.stringency} --basename {params.basename} -j 7 --paired --gzip -o {params.output_dir} {input} > {log} 2>&1"


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
		Lmax=0.04
	message: "RNA ALIGN: Align {wildcards.sample} to the genome by STAR"
	shell:
		"""STAR --runThreadN {params.thread} --runMode alignReads --genomeDir {params.ref} --readFilesCommand zcat --readFilesIn {input[0]} {input[1]} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
		--outFilterMismatchNmax {params.Nmax} --outFilterMismatchNoverLmax {params.Lmax} """



rule feature_count:
	input: 	
		lambda wildcards: expand(config['rnadir'] + config['bam'] + "{sample}Aligned.sortedByCoord.out.bam", sample = get_dnausecells(wildcards))
	output: 
		temp(config['rnadir'] + config['bam'] + "count.txt"),
		temp(config['rnadir'] + config['bam'] + "count.txt.summary"),
		config['rnadir'] + config['bam'] + "rawcount.txt",
	params:
		featuretype="exon",
		ref_gtf=config['star_annotation'],
		threads=5,
		attributetype="gene_id"
	message: "QUANTIFICATION: RNA quantification using FeatureCount"
	shell:
		"""featureCounts -T {params.threads} -p -t {params.featuretype} -g {params.attributetype} -a {params.ref_gtf} -o {output[0]} {input} &&\
		cut -f 1,7- {output[0]} > {output[2]} &&\
		sed -i '1d' {output[2]}"""








