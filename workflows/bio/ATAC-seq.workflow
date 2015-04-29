# -*- snakemake -*-
import os
import pysam
from snakemakelib.io import set_output
from snakemakelib.config import update_sml_config, get_sml_config
from snakemakelib.bio.ngs.targets import generic_target_generator

def _merge_suffix():
    """Determine the merge suffix of the run files"""
    atac_cfg = get_sml_config('workflows.bio.atac_seq')
    if atac_cfg['trimadaptor']:
        return ".trimmed.sort.bam"
    else:
        return ".sort.bam"

# Default configuration settings custom-tailored for ATAC-Seq analysis
atac_config = {
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie',
        'peakcallers' : ['zinba', 'dfilter', 'macs2'],
        'trimadaptor' : True,
    },
    'settings' : {
        'temp_rules' : [],
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : '.sort.bam',
        },
    },
    'bio.ngs.enrichment.macs' : {
        'callpeak' : {
            'options' : '-g hs --nomodel --shift 37 --extsize 73 -q 0.01',
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

ALIGN_TARGET_SUFFIX = ".bam"
if atac_config['workflows.bio.atac_seq']['trimadaptor']:
    atac_config['bio.ngs.qc.picard']['merge_sam']['suffix'] = '.trimmed.sort.bam'
    ALIGN_TARGET_SUFFIX = ".trimmed.bam"

aligner = atac_config['workflows.bio.atac_seq']['aligner']
key = 'bio.ngs.align.' + aligner

update_sml_config(atac_config)
update_sml_config({key : aligner_config[key]})

ngs_cfg = get_sml_config('bio.ngs.settings')
main_cfg = get_sml_config('settings')
atac_cfg = get_sml_config('workflows.bio.atac_seq')

p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: os.path.join(p, "bio/ngs/align", "blat.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(p, "bio/ngs/enrichment", "zinba.rules")
include: os.path.join(p, "bio/ngs/enrichment", "dfilter.rules")
include: os.path.join(p, "bio/ngs/enrichment", "macs.rules")
if atac_cfg['trimadaptor']:
    include: os.path.join(p, "bio/ngs/qc", "cutadapt.rules")

ruleorder: picard_merge_sam > picard_sort_bam 
ruleorder: picard_sort_bam > picard_add_or_replace_read_groups
ruleorder: picard_add_or_replace_read_groups > picard_mark_duplicates
ruleorder: picard_mark_duplicates > atacseq_correct_coordinates
ruleorder: atacseq_correct_coordinates > bowtie_align
ruleorder: picard_sort_bam > bowtie_align
ruleorder: picard_merge_sam > bowtie_align
ruleorder: picard_mark_duplicates > bowtie_align
ruleorder: dfilter_run_dfilter_bam > bedtool_bamtobed
ruleorder: macs_callpeak_treatment_only_bam > bedtool_bamtobed

# Set temporary and protected outputs
set_output(workflow,
           temp_rules = main_cfg['temp_rules'] + main_cfg['temp_rules_default'],
           temp_filetypes=main_cfg['temp_filetypes'] + main_cfg['temp_filetypes_default'],
           protected_rules = main_cfg['protected_rules'] + main_cfg['protected_rules_default'],
           protected_filetypes=main_cfg['protected_filetypes'] + main_cfg['protected_filetypes_default'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

ALIGN_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].run_id_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = ALIGN_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = MERGE_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg)

ZINBA_TARGET_SUFFIX = ".sort.merge.offset.zinba.peaks"
ZINBA_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = ZINBA_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg) if 'zinba' in atac_cfg['peakcallers'] else []

DFILTER_TARGET_SUFFIX = ".sort.merge.offset.dfilt.bed"
DFILTER_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = DFILTER_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg) if 'dfilter' in atac_cfg['peakcallers'] else []

MACS2_TARGET_SUFFIX = ".sort.merge.offset_peaks.bed"
MACS2_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, src_re = ngs_cfg['sampleorg'].raw_run_re, target_suffix = MACS2_TARGET_SUFFIX, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'],  **ngs_cfg) if 'macs2' in atac_cfg['peakcallers'] else []

DUP_METRICS_SUFFIX=".sort.merge.dup.dup_metrics"
DUP_METRICS_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, target_suffix =  DUP_METRICS_SUFFIX, src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], **ngs_cfg)

ALIGN_METRICS_SUFFIX=".sort.merge.dup.align_metrics"
ALIGN_METRICS_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, target_suffix =  ALIGN_METRICS_SUFFIX, src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], **ngs_cfg)

INSERT_METRICS_SUFFIX=".sort.merge.dup.insert_metrics"
INSERT_METRICS_TARGETS = generic_target_generator(tgt_re = ngs_cfg['sampleorg'].sample_re, target_suffix =  INSERT_METRICS_SUFFIX, src_re = ngs_cfg['sampleorg'].raw_run_re, filter_suffix = ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix'], **ngs_cfg)


# Rules
rule atacseq_align:
    """Run ATAC-seq alignment"""
    input: ALIGN_TARGETS

rule atacseq_merge:
    """Run ATAC-seq alignment, duplication removal and merge"""
    input: MERGE_TARGETS

rule atacseq_metrics:
    """Run ATAC-seq alignment and corresponding metrics only"""
    input: DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS

rule atacseq_all:
    """Run ATAC-seq pipeline"""
    input: DFILTER_TARGETS + ZINBA_TARGETS + MACS2_TARGETS + DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS

rule atacseq_correct_coordinates:
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
            # Modify s here
            if not s.is_unmapped:
                l = samfile.lengths[s.rname]
                if not s.is_reverse:
                    s.pos = min(l, s.pos + 4)
                    s.pnext = max(0, s.pnext - 5)
                else:
                    s.pos = max(0, s.pos - 5)
                    s.pnext = min(l, s.pnext + 4)
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

# Other papers
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004994
# use macs2 with -g dm –nomodel –shiftsize 50 –q 0.01; however single-end 50bp reads

