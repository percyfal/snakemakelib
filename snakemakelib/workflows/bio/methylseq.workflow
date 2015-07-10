# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.utils import rreplace
from snakemakelib.bio.ngs.methylseq.bismark import report_label, align_suffix
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.regexp import SampleRegexp

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
    """
    ngs_cfg = get_sml_config('bio.ngs.settings')
    picard_cfg = get_sml_config('bio.ngs.qc.picard')
    sources = generic_target_generator(tgt_re=ngs_cfg['sampleorg'].run_id_re, target_suffix = align_suffix(), src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'] + "$", **ngs_cfg)
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

FASTQC_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = "_1_fastqc/fastqc_report.html", src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'] + "$", **ngs_cfg) + \
    generic_target_generator(tgt_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = "_2_fastqc/fastqc_report.html", src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'] + "$", **ngs_cfg)

sr = SampleRegexp (os.path.join(os.path.dirname(ngs_cfg['sampleorg'].sample_re.pattern), "CpG_OB_" + os.path.basename(ngs_cfg['sampleorg'].sample_re.pattern)))
BISMARK_TARGETS = generic_target_generator(tgt_re = sr, target_suffix = ".merge.deduplicated.txt.gz", src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], **ngs_cfg)

BISMARK_REPORT_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, target_suffix = ".merge.deduplicated.bam{report_label}.html".format(report_label=report_label()), src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], **ngs_cfg)

# All rules
rule bismark_all:
    """Run all the analyses"""
    input: FASTQC_TARGETS + BISMARK_TARGETS + BISMARK_REPORT_TARGETS

rule run_bismark_fastqc:
    """Fastqc target rule. Run fastqc on files defined in FASTQC_TARGETS"""
    input: FASTQC_TARGETS

rule run_bismark_run:
    """bismark target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_TARGETS

rule run_bismark_report:
    """bismark report target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_REPORT_TARGETS

rule bismark_targets:
    """List currently defined targets"""
    run:
      print ("Fastqc targets: ", FASTQC_TARGETS)
      print ("Bismark targets: ", BISMARK_TARGETS)
      print ("Bismark report targets: ", BISMARK_REPORT_TARGETS)
