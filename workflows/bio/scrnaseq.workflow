# -*- snakemake -*-
import os
from snakemake.workflow import workflow
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.utils import ReadGroup, find_files

def _merge_suffix(aligner, quantification=[]):
    align_cfg = get_sml_config('bio.ngs.align.' + aligner)
    if aligner == "star":
        if "rsem" in quantification:
            return ".Aligned.toTranscriptome.out_unique.bam"
        return align_cfg['align']['suffix'].replace('.bam', '_unique.bam')
    elif aligner in ["bowtie", "bowtie2"]:
        return '_unique.bam'

def find_scrnaseq_merge_inputs(wildcards):
    """Find platform unit specific aligned bam files as input to picard
merge.
    """
    ngs_cfg = get_sml_config('bio.ngs.settings')
    picard_cfg = get_sml_config('bio.ngs.qc.picard')
    rg = ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], cfg=ngs_cfg, path=wildcards.path)
    fmt = ngs_cfg['run_id_pfx_fmt'] + _merge_suffix(ngs_cfg['aligner'], ngs_cfg['rnaseq']['quantification'])
    sources = generic_target_generator(fmt=os.path.basename(fmt), rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=wildcards.path)
    sources = [src for src in sources if os.path.basename(src).startswith(wildcards.prefix)]
    return sources

config_default = {
    'bio.ngs.settings' : {
        'aligner' : 'bowtie',
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'label' : 'merge',
            'inputfun' : find_scrnaseq_merge_inputs,
        },
    },
}


update_sml_config(config_default)
ngs_cfg = get_sml_config('bio.ngs.settings')
aligner = ngs_cfg['aligner']

# Include necessary snakemakelib rules
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
align_include = os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: align_include
include: os.path.join(p, "bio/ngs/tools", "bamtools.rules")
include: os.path.join(p, "bio/ngs/qc", "rseqc.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/tools", "samtools.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "rpkmforgenes.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "rsem.rules")
if aligner in ["bowtie", "bowtie2"]:
    ruleorder: bamtools_filter > picard_merge_sam > picard_sort_bam > bowtie_align

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'scrnaseq.workflow'; set workdir before include statement!")
path = workflow._workdir

ALIGN_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + _merge_suffix(aligner), rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

RSEQC_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge_rseqc/rseqc_qc_8.txt', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

RPKMFORGENES_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge.rpkmforgenes', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path) if 'rpkmforgenes' in ngs_cfg['rnaseq']['quantification']  else []

RSEM_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge.sort.transcript.bam', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path) if 'rsem' in ngs_cfg['rnaseq']['quantification']  else []

# All rules
rule scrnaseq_all:
    """Run scRNAseq pipeline"""
    input: ALIGN_TARGETS + RSEQC_TARGETS + RPKMFORGENES_TARGETS + RSEM_TARGETS

rule scrnaseq_align:
    """Run alignments"""
    input: ALIGN_TARGETS

rule scrnaseq_rseqc:
    """Run RSeQC"""
    input: RSEQC_TARGETS

rule scrnaseq_rpkmforgenes:
    """Run rpkmforgenes"""
    input: RPKMFORGENES_TARGETS

rule scrnaseq_rsem:
    """Run rpkmforgenes"""
    input: RSEM_TARGETS
    
rule scrnaseq_targets:
    """Print targets """
    run:
        print (ALIGN_TARGETS)
        print (RSEQC_TARGETS)
        print (RPKMFORGENES_TARGETS)
        print (RSEM_TARGETS)

rule clean:
    """Clean working directory. WARNING: will remove all files except
    (.fastq|.fastq.gz) and csv files
    """
    params: d = workflow._workdir
    shell: 'for f in `find  {params.d} -type f -name "*" | grep -v ".fastq$" | grep -v ".fastq.gz$" | grep -v ".csv$"`; do echo removing $f; rm -f $f; done'
