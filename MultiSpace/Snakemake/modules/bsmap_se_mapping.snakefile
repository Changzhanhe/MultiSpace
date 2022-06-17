def bsmap_target(wildcards):
	ls = []
	for sample in config['runs']:
		ls.append(config['dnadir'] + config['bam'] + "%s/%s.sort.rmdup.bam" % (sample,sample))
		ls.append(config['dnadir'] + config['bam'] + "%s/%s.sort.rmdup.bam.bai" % (sample,sample))
	return ls


def getAlignFastq(wildcards):
	tmp = []
	tmp.append(config['dnadir'] + config['fastq'] + "%s/%s_val_1.fq.gz" % (wildcards.sample,wildcards.sample))
	tmp.append(config['dnadir'] + config['fastq'] + "%s/%s_val_2.fq.gz" % (wildcards.sample,wildcards.sample))
	return tmp


rule bsmap_mapping:
	input: getAlignFastq
	output:
		config['dnadir'] + config['bam'] + "{sample}/{sample}_1.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}_2.bam"
	params:
		ref_fa=config['lambda_fasta'],
		map_str=1
	message: "COOL-seq ALIGN: Align {wildcards.sample} to the genome by bsmap"
	log: config['dnadir'] + config['bam'] + "{sample}/{sample}.log"
	shell:
		"bsmap -p 10 -a {input[0]} -d {params.ref_fa} -n {params.map_str}  -o {output[0]} > {log} 2>&1 &&\
		 bsmap -p 10 -a {input[1]} -d {params.ref_fa} -n {params.map_str}  -o {output[1]} >> {log} 2>&1"

rule bam_rmdup:
	input:
		config['dnadir'] + config['bam'] + "{sample}/{sample}_1.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}_2.bam",
	output:
		config['dnadir'] + config['bam'] + "{sample}/{sample}.1.sort.rmdup.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}.2.sort.rmdup.bam",
	params:
		maxmem=1000000000
	message: "REMOVE DUP: remove {wildcards.sample} bam duplicate using samtools"
	shell:
		"samtools view -@ 8 -u {input[0]} | samtools sort -m {params.maxmem} - | samtools rmdup -s - {output[0]} && samtools view -@ 8 -u {input[1]} |\
		 samtools sort -m {params.maxmem} - | samtools rmdup -s - {output[1]}" 


rule merge_2bam:
	input:
		config['dnadir'] + config['bam'] + "{sample}/{sample}.1.sort.rmdup.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}.2.sort.rmdup.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}_1.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}_2.bam",
	output:
		temp(config['dnadir'] + config['bam'] + "{sample}/{sample}.rmdup.bam"),
		config['dnadir'] + config['bam'] + "{sample}/{sample}.sort.rmdup.bam",
		config['dnadir'] + config['bam'] + "{sample}/{sample}.sort.rmdup.bam.bai",		
	params:
		maxmem=1000000000
	message: "MERGE 3 BAM: merge sort and get {wildcards.sample} bam files index"
	shell:
		"samtools merge -f {output[0]} {input[0]} {input[1]} && samtools sort -@ 8 -m {params.maxmem} {output[0]} -o {output[1]} && samtools index {output[1]} \
		&& rm {input[0]} {input[1]} {input[2]} {input[3]}"


