# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.utils import ReadGroup

atac_config = {
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie',
    },
    'bio.ngs.settings' : {
        'sample_organization' : 'Illumina@SciLife',
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


ruleorder: picard_add_or_replace_read_groups > samtools_index

ngs_cfg = get_sml_config('bio.ngs.settings')
if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")
path = workflow._workdir

# Target suffixes
BAM_TARGET_SUFFIX = ".sort.merge.rg.dup.bam"

BAM_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + BAM_TARGET_SUFFIX, rg = ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

# Rules
rule atacseq_finaldup:
    """Run ATAC-seq alignment, merge and duplication removal"""
    input: BAM_TARGETS
