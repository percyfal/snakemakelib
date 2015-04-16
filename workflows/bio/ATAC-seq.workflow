# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator

atac_config = {
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie',
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

ruleorder: picard_add_or_replace_read_groups > bowtie_align > samtools_index > picard_build_bam_index

ruleorder: picard_mark_duplicates > bowtie_align

ngs_cfg = get_sml_config('bio.ngs.settings')

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

# Target suffixes
ALIGN_TARGET_SUFFIX = ".sort.merge.rg.dup.bam"
ALIGN_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].run_id_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = ALIGN_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

# Rules
rule atacseq_finaldup:
    """Run ATAC-seq alignment, merge and duplication removal"""
    input: ALIGN_TARGETS

