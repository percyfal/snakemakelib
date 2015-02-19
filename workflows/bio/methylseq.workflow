# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.methylseq.bismark import report_label, align_suffix

def find_report_inputs(wildcards):
    """Find bismark align report files to use as input to
    bismark_merge_alignment_reports.

    Assumes folder structure:

    {path}/{flowcell}/{prefix}_[PE|SE]_report.txt
    """
    subdirs = [d for d in os.listdir(wildcards.path) if os.path.isdir(os.path.join(wildcards.path, d))]
    sources = []
    for d in subdirs:
        reportsrc = sorted(glob.glob("{path}{sep}*{suffix}".format(path=os.path.join(wildcards.path, d), sep=os.sep, suffix=report_label() + ".txt")))
        sources += reportsrc
    return sources

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
        bismarkbam = [x.replace(sfx, "") + align_suffix() for x in fqsrc if x.endswith(sfx)]
        sources += bismarkbam
    return sources

methylation_config = {
    'bio.ngs.methylseq.bismark' : {
        'methXtract' : {
            'options' : "--ignore_r2 2 --counts  --gzip -p --no_overlap",
        },
        'report' : {
            'inputfun' : find_report_inputs,
        },
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : align_suffix,
            'label' : 'merge',
            'inputfun' : find_meth_merge_inputs,
            'options' : "SORT_ORDER=queryname", # methylation extraction requires queries to be adjacent
        },
    },
}

update_sml_config(methylation_config)

# Include necessary snakemakelib rules
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(p, "bio/ngs/methylseq", "bismark.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs", "settings.rules")

# Get relevant config sections
bismark_cfg = get_sml_config('bio.ngs.methylseq.bismark')
qc_cfg = get_sml_config('bio.ngs.qc.sequenceprocessing')
cfg = get_sml_config('bio.ngs.settings')
path = cfg.get('path') if not cfg.get('path') is None else os.curdir

FASTQC_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_1_fastqc.html {path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_2_fastqc.html".split(), sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

BISMARK_TARGETS = expand("{path}/{sample}/CpG_OB_{sample}.merge.deduplicated.txt.gz", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

BISMARK_REPORT_TARGETS = expand("{path}/{sample}/{sample}.merge.deduplicated.bam{report_label}.html", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path, report_label=report_label())

# All rules
rule all:
    """Run all the analyses"""
    input: FASTQC_TARGETS + BISMARK_TARGETS + BISMARK_REPORT_TARGETS

rule run_fastqc:
    """Fastqc target rule. Run fastqc on files defined in FASTQC_TARGETS"""
    input: FASTQC_TARGETS

rule run_bismark:
    """bismark target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_TARGETS

rule run_bismark_report:
    """bismark report target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_REPORT_TARGETS

rule list_targets:
    """List currently defined targets"""
    run:
      print ("Sample targets: ", SAMPLE_TARGETS)
      print ("Fastqc targets: ", FASTQC_TARGETS)
      print ("Bismark targets: ", BISMARK_TARGETS)
      print ("Bismark report targets: ", BISMARK_REPORT_TARGETS)
