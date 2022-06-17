def trim_target(wildcards):
    ls = []
    for sample in config['runs']:
        ls.append(config['dnadir'] + config['fastq'] + "%s/%s_val_1.fq.gz" % (sample,sample))
        ls.append(config['dnadir'] + config['fastq'] + "%s/%s_val_2.fq.gz" % (sample,sample))
    return ls


def get_fastq(wildcards):
    return glob.glob(config['dnadir'] + config['raw'] + wildcards.sample + "_[1-2].fastq.gz")


rule trim_adapter:
    input: get_fastq
    output:
        config['dnadir'] + config['fastq'] + "{sample}/{sample}_val_1.fq.gz",
        config['dnadir'] + config['fastq'] + "{sample}/{sample}_val_2.fq.gz",
        temp(config['dnadir'] + config['fastq'] + "{sample}/{sample}_1.fastq.gz_trimming_report.txt"),
        temp(config['dnadir'] + config['fastq'] + "{sample}/{sample}_2.fastq.gz_trimming_report.txt")
    params:
        basename=lambda wildcards: "%s" % (wildcards.sample),
        quality=20,
        stringency=3,
        length=50,
        clip_R1=6,
        clip_R2=6,
        output_dir=lambda wildcards: config['dnadir'] + config['fastq'] + "%s/" % (wildcards.sample)
    message: "TRIM: Trim adaptors for paired-end using trim_galore - {wildcards.sample}"
    log: config['dnadir'] + config['fastq'] +  "log/{sample}.log"
    shell: 
        "trim_galore --quality {params.quality} --stringency {params.stringency} --length {params.length} \
        --clip_R1 {params.clip_R1} --basename {params.basename} --clip_R2 {params.clip_R2} -j 7 --paired --trim1 --phred33 --gzip -o {params.output_dir} {input} > {log} 2>&1"



