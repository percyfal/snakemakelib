# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, sml_rules_path, get_sml_config

def find_meth_merge_inputs(wildcards):
    """Find bismark aligned bam files as input to picard merge.

    Assumes folder structure:

    {path}/{flowcell}/{prefix}.bam

    to be merged into {path}, which often represents sample.
    """
    subdirs = [d for d in os.listdir(wildcards.path) if os.path.isdir(os.path.join(wildcards.path, d))]
    sources = []
    for d in subdirs:
        fqsrc = sorted(glob.glob("{path}{sep}*{suffix}".format(path=os.path.join(wildcards.path, d), sep=os.sep, suffix=sml_config['bio.ngs.settings']['fastq_suffix'])))
        sfx = sml_config['bio.ngs.settings']['read1_label'] + sml_config['bio.ngs.settings']['fastq_suffix']
        bismarkbam = [x.replace(sfx, "") + sml_config['bio.ngs.methylseq.bismark']['align']['suffix'] for x in fqsrc if x.endswith(sfx)]
        sources += bismarkbam
    return sources

methylation_config = {
    'bio.ngs.methylseq.bismark' : {
        'methXtract' : {
            'options' : "--ignore_r2 2 --counts  --gzip -p --no_overlap",
        },
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : '_bismark_bt2_pe.bam',
            'label' : 'merge',
            'inputfun' : find_meth_merge_inputs,
            'options' : "SORT_ORDER=queryname", # methylation extraction requires queries to be adjacent
        },
    },
}

update_sml_config(methylation_config)

include: os.path.join(sml_rules_path(), 'settings.rules')
include: os.path.join(sml_rules_path(), 'utils.rules')
include: os.path.join(sml_rules_path(), "bio/ngs", "settings.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/methylseq", "bismark.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/qc", "picard.rules")

bismark_cfg = get_sml_config('bio.ngs.methylseq.bismark')
qc_cfg = get_sml_config('bio.ngs.qc.sequenceprocessing')
cfg = get_sml_config('bio.ngs.settings')

path = cfg.get('path') if not cfg.get('path') is None else os.curdir

FASTQC_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_1_fastqc.html {path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_2_fastqc.html".split(), sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

BISMARK_TARGETS = expand("{path}/{sample}/CpG_OB_{sample}.merge.deduplicated.txt.gz", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

#BISMARK_REPORT_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_bismark_bt2_PE_report.html", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

#BISMARK_REPORT_TARGETS = expand("{path}/{sample}/{sample}_PE_report.html", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

# All rules
rule all:
    """Run all the analyses"""
    input: FASTQC_TARGETS + BISMARK_TARGETS# + BISMARK_REPORT_TARGETS

rule run_fastqc:
    """Fastqc target rule. Run fastqc on files defined in FASTQC_TARGETS"""
    input: FASTQC_TARGETS

rule run_bismark:
    """bismark target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_TARGETS

# TODO: generic rule
rule list_targets:
    """List currently defined targets"""
    run:
      print ("In list_targets")
      print (SAMPLE_TARGETS)
      print (FASTQC_TARGETS)
      print (BISMARK_TARGETS)
      print (BISMARK_REPORT_TARGETS)

# Redefine this rule for now as the input are files from trim_galore
rule methylseq_bismark_align:
    """methylseq.workflow: Run bismark align."""
    params: options = bismark_cfg['align']['options'],
            cmd = bismark_cfg['align']['cmd'],
            ref = bismark_cfg['ref'],
    threads: bismark_cfg['align']['threads']
    input: "{path}" + os.sep + "{prefix}" + ngs_cfg['read1_label'] + qc_cfg['trim_galore']['read1_suffix'],\
           "{path}" + os.sep + "{prefix}" + ngs_cfg['read2_label'] + qc_cfg['trim_galore']['read2_suffix']
    output: "{path}" + os.sep + "{prefix}" + bismark_cfg['align']['suffix'], "{path}" + os.sep + "{prefix}_PE_report.txt"
    shell: "{params.cmd} {params.options} {params.ref} -1 {input[0]} -2 {input[1]} -o {wildcards.path} -p {threads} -B {wildcards.prefix}"
