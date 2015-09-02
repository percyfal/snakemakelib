# -*- snakemake -*-
include: "../settings.rules"

config_default = {
    'bio.ngs.qc.trim_galore' : {
        'cmd' : "trim_galore",
        'options' : "--paired --trim1 --phred33 --fastqc",
        'read1_suffix' : "_val_1.fq.gz",
        'read2_suffix' : "_val_2.fq.gz",
    },
}

update_config(config_default, config)
config = config_default


rule trim_galore_main:
    """trim_galore: run trim_galore on a fastq file"""
    params: cmd = config['bio.ngs.qc.trim_galore']['cmd'],
            options = " ".join([config['bio.ngs.qc.trim_galore']['options'],"--gzip"])
    input: "{prefix}" + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix'],\
           "{prefix}" + config['bio.ngs.settings']['read2_label'] + config['bio.ngs.settings']['fastq_suffix']
    output: "{prefix}" + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.qc.trim_galore']['read1_suffix'],\
            "{prefix}" + config['bio.ngs.settings'] ['read2_label'] + config['bio.ngs.qc.trim_galore']['read2_suffix']
    shell: "{params.cmd} {params.options} {input} -o $(dirname {wildcards.prefix})"
