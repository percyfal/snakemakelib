# -*- snakemake -*-
import os
import pysam
from snakemakelib.io import set_temp_output
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator

atac_config = {
    # 'settings' : {
    #     'temp_rules' : ['sratools_prefetch', 'sratools_fastq_dump']
    # },
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
        'version2' : False,
        'bowtie' : {
            'options' : '-X 2000 -m1 --chunkmbs 1000',
        },
    },
    'bio.ngs.align.bwa' : {
    },
}
aligner = atac_config['workflows.bio.atac_seq']['aligner']
key = 'bio.ngs.align.' + aligner

update_sml_config(atac_config)
update_sml_config({key : aligner_config[key]})

p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: os.path.join(p, "bio/ngs/align", "blat.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/enrichment", "zinba.rules")

ruleorder: picard_merge_sam > picard_sort_bam > picard_add_or_replace_read_groups > picard_mark_duplicates > atacseq_correct_coordinates_for_zinba > bowtie_align

#ruleorder: ucsc_autosome_reference > peakseq_mappability_write_chromosome

ngs_cfg = get_sml_config('bio.ngs.settings')
main_cfg = get_sml_config('settings')

# Set temporary outputs
set_temp_output(workflow, main_cfg['temp_rules'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = MERGE_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

ZINBA_TARGET_SUFFIX = ".sort.merge.offset.zinba"
ZINBA_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = ZINBA_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

# Rules
rule atacseq_merge:
    """Run ATAC-seq alignment, duplication removal and merge"""
    input: MERGE_TARGETS

rule atacseq_all:
    """Run ATAC-seq pipeline"""
    input: ZINBA_TARGETS + ['zinba_alignability/']

rule atacseq_correct_coordinates_for_zinba:
    """From Buenrostro paper: 

    'Previous descriptions of the Tn5 transposase show that the
    transposon binds as a dimer and inserts two adaptors separated by
    9 bp (ref. 11). Therefore, all reads aligning to the + strand were
    offset by +4 bp, and all reads aligning to the – strand were
    offset −5 bp'
    """
    input: bam = "{prefix}.bam"
    output: bam = "{prefix}.offset.bam"
    run:
        # Use pysam to modify input
        samfile = pysam.AlignmentFile(input.bam, "rb")
        outfile = pysam.AlignmentFile(output.bam, "wb", template=samfile)
        for s in samfile:
            # Modify s here NB: currently not checking if position is
            # <0 or >chr_len, should be done
            if not s.is_unmapped:
                if not s.is_reverse:
                    s.pos += 4
                    s.pnext -= 5
                else:
                    s.pos -= 5
                    s.pnext += 4
            outfile.write(s)


