# -*- snakemake -*-
import os
import pysam
from snakemakelib.io import set_temp_output
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator

atac_config = {
    'settings' : {
        'temp_rules' : [],
        'temp_rules_default' : ['sratools_prefetch'],
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
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: os.path.join(p, "bio/ngs/align", "blat.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/enrichment", "zinba.rules")
include: os.path.join(p, "bio/ngs/enrichment", "dfilter.rules")

ruleorder: picard_merge_sam > picard_sort_bam 
ruleorder: picard_sort_bam > picard_add_or_replace_read_groups
ruleorder: picard_add_or_replace_read_groups > picard_mark_duplicates
ruleorder: picard_mark_duplicates > atacseq_correct_coordinates_for_zinba
ruleorder: atacseq_correct_coordinates_for_zinba > bowtie_align
ruleorder: picard_sort_bam > bowtie_align
ruleorder: picard_merge_sam > bowtie_align
ruleorder: picard_mark_duplicates > bowtie_align
ruleorder: dfilter_run_dfilter_bam > bedtool_bamtobed

ngs_cfg = get_sml_config('bio.ngs.settings')
main_cfg = get_sml_config('settings')

# Set temporary outputs
set_temp_output(workflow, main_cfg['temp_rules'] + main_cfg['temp_rules_default'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = MERGE_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

ZINBA_TARGET_SUFFIX = ".sort.merge.offset.zinba.peaks"
ZINBA_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = ZINBA_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

DFILTER_TARGET_SUFFIX = ".sort.merge.offset.dfilt.bed"
DFILTER_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = DFILTER_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

# Rules
rule atacseq_merge:
    """Run ATAC-seq alignment, duplication removal and merge"""
    input: MERGE_TARGETS

rule atacseq_all:
    """Run ATAC-seq pipeline"""
    input: DFILTER_TARGETS + ZINBA_TARGETS

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

#
# Putative additional data and methods
#
# Section ATAC-seq insertion size enrichment analysis
# No Segmentation data available for current ensembl build ?!? Use old data:
# ftp://ftp.ensembl.org/pub/release-73/regulation/homo_sapiens/Segmentation_GM12878.gff.gz
#
# Section Nucleosome positioning
# Danpos, Dantools: https://sites.google.com/site/danposdoc/
# 
# Section ChIP-seq peak-calling and clustering
# ChIP data for 50 antibodies downloaded from
# Stanford/Yale/USC/Harvard (SYDH) ENCODE data repository available at the UCSC genome browser:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/
# 
# Section Footprinting using CENTIPEDE

# Use motif data from http://compbio.mit.edu/encode-motifs/; most
# likely need matches.txt.gz to find genomic regions matching a motif;
# the genome-wide set of motifs is in motifs.txt (?).
