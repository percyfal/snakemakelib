# -*- snakemake -*-
import os
from snakemake.workflow import workflow
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.utils import ReadGroup, find_files

def find_scrnaseq_merge_inputs(wildcards):
    """Find platform unit specific aligned bam files as input to picard
merge.
    """
    cfg = get_sml_config('bio.ngs.settings')
    picard_cfg = get_sml_config('bio.ngs.qc.picard')
    rg = ReadGroup(cfg['run_id_pfx_re'] + cfg['read1_label'] + cfg['fastq_suffix'], cfg=cfg, path=wildcards.path)
    fmt = cfg['run_id_pfx_fmt'] + picard_cfg['merge_sam']['suffix']
    inputs = find_files(path=wildcards.path, re_str=rg.pattern)
    if inputs:
        # FIXME: Yet another problem caused by readgroup parsing
        # inconsistency. PATH should always be included in read group
        # information

        #sources = [fmt.format(PATH=workflow._workdir, **dict(rg.parse(f))) for f in inputs]
        sources = [fmt.format(**dict(rg.parse(f))) for f in inputs]
    return sources

config_default = {
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : ".Aligned.out_unique.bam",
            'label' : 'merge',
            'inputfun' : find_scrnaseq_merge_inputs,
            'options' : "SORT_ORDER=coordinate", 
        },
    },
}

update_sml_config(config_default)    

# Include necessary snakemakelib rules
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs/align", "star.rules")
include: os.path.join(p, "bio/ngs/qc", "rseqc.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/tools", "bamtools.rules")
include: os.path.join(p, "bio/ngs/tools", "samtools.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "rpkmforgenes.rules")

ngs_cfg = get_sml_config('bio.ngs.settings')
if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'scrnaseq.workflow'; set workdir before include statement!")
path = workflow._workdir

STAR_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.Aligned.out_unique.bam', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

#RSEQC_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.Aligned.out_unique_rseqc/rseqc_qc_8.txt', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)
RSEQC_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge.sort_rseqc/rseqc_qc_8.txt', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

#RPKMFORGENES_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.Aligned.out_unique.rpkmforgenes', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)
RPKMFORGENES_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge.sort.rpkmforgenes', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

# All rules
rule scrnaseq_all:
    """Run scRNAseq pipeline"""
    input: STAR_TARGETS + RSEQC_TARGETS + RPKMFORGENES_TARGETS

rule scrnaseq_star:
    """Run STAR alignments"""
    input: STAR_TARGETS

rule scrnaseq_rseqc:
    """Run RSeQC"""
    input: RSEQC_TARGETS

rule scrnaseq_rpkmforgenes:
    """Run rpkmforgenes"""
    input: RPKMFORGENES_TARGETS
    
rule scrnaseq_targets:
    """Print targets """
    run:
        print (STAR_TARGETS)
        print (RSEQC_TARGETS)
        print (RPKMFORGENES_TARGETS)

rule clean:
    """Clean working directory. WARNING: will remove all files except
    (.fastq|.fastq.gz) and csv files
    """
    params: d = workflow._workdir
    shell: 'for f in `find  {params.d} -type f -name "*" | grep -v ".fastq$" | grep -v ".fastq.gz$" | grep -v ".csv$"`; do echo removing $f; rm -f $f; done'
