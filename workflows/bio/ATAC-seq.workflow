# -*- snakemake -*-
import os
from snakemakelib.io import set_temp_output
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator

atac_config = {
    'settings' : {
        'temp_rules' : ['sratools_prefetch', 'sratools_fastq_dump']
    },
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie',
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : '.sort.dup.bam',
        },
    },
}
aligner_config = {
    'bio.ngs.align.bowtie' : {
        'bowtie' : {
            'options' : '-X 2000 -m1',
        },
    },
    'bio.ngs.align.bwa' : {
    },
}
aligner = atac_config['workflows.bio.atac_seq']['aligner']
key = 'bio.ngs.align.' + aligner

update_sml_config(atac_config)
update_sml_config(aligner_config[key])

p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/enrichment", "zinba.rules")

ruleorder: picard_merge_sam > picard_sort_bam > picard_add_or_replace_read_groups > picard_mark_duplicates > bowtie_align

ngs_cfg = get_sml_config('bio.ngs.settings')
main_cfg = get_sml_config('settings')

# Set temporary outputs
set_temp_output(workflow, main_cfg['temp_rules'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = MERGE_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

# Rules
rule atacseq_finaldup:
    """Run ATAC-seq alignment, merge and duplication removal"""
    input: MERGE_TARGETS

rule atacseq_all:
    """Run ATAC-seq pipeline"""
    input: MERGE_TARGETS

