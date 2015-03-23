# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.utils import rreplace
from snakemakelib.bio.ngs.methylseq.bismark import report_label, align_suffix
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.utils import ReadGroup

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
            'options' : "--ignore_r2 2 --counts  --gzip -p --no_overlap --bedGraph --report",
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

FASTQC_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + "_1_fastqc.html", rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

BISMARK_TARGETS = generic_target_generator(fmt=rreplace(ngs_cfg['sample_pfx_fmt'], "{SM}", "CpG_OB_{SM}", 1) + ".merge.deduplicated.txt.gz", rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

BISMARK_REPORT_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + ".merge.deduplicated.bam{report_label}.html".format(report_label=report_label()), rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

# All rules
rule bismark_all:
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
      print ("Fastqc targets: ", FASTQC_TARGETS)
      print ("Bismark targets: ", BISMARK_TARGETS)
      print ("Bismark report targets: ", BISMARK_REPORT_TARGETS)
