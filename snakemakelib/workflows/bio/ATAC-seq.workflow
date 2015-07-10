# -*- snakemake -*-
import os
import pysam
from jinja2 import Environment, PackageLoader
from snakemakelib.io import set_output
from snakemakelib.config import update_snakemake_config
from snakemakelib.bio.ngs.targets import generic_target_generator

##############################
# Functions
##############################
def _merge_suffix():
    """Determine the merge suffix of the run files"""
    atac_cfg = config['workflows.bio.atac_seq']
    if atac_cfg['trimadaptor']:
        return ".trimmed.sort.bam"
    else:
        return ".sort.bam"


##############################
# Default configuration settings custom-tailored for ATAC-Seq analysis
##############################
atac_config = {
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie',
        'peakcallers' : ['dfilter', 'macs2'],
        'trimadaptor' : True,
        'bamfilter' : True,
    },
    'settings' : {
        'temp_rules' : [],
    },
    'bio.ngs.qc.picard' : {
        'merge_sam' : {
            'suffix' : '.sort.bam',
        },
        'qcrules' : ['picard_collect_insert_size_metrics',
                     'picard_mark_duplicates'],
    },
    'bio.ngs.enrichment.macs' : {
        'callpeak' : {
            'options' : '-g hs --nomodel --nolambda --keep-dup all --call-summits -B',
        },
    },
    'bio.ngs.tools.bamtools' : {
        'filter' : {
            'options' : {'mapQuality': ">=30",
                         'isProperPair': 'true'},
        },
    },
}
aligner_config = {
    'bio.ngs.align.bowtie' : {
        'version2' : True,
        'bowtie' : {
            'options' : '-X 2000',
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

config = update_snakemake_config(config, atac_config)
config = update_snakemake_config(config, {key :  aligner_config[key]})
ngs_cfg = config['bio.ngs.settings']
main_cfg = config['settings']
atac_cfg = config['workflows.bio.atac_seq']

##############################
# Include statements
##############################
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/align", aligner + ".rules")
include: os.path.join(p, "bio/ngs/align", "blat.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(p, "bio/ngs/qc", "qualimap.rules")
include: os.path.join(p, "bio/ngs/chromatin", "danpos.rules")
if 'dfilter' in config['workflows.bio.atac_seq']['peakcallers']:
    include: os.path.join(p, "bio/ngs/enrichment", "dfilter.rules")
if 'macs2' in config['workflows.bio.atac_seq']['peakcallers']:
    include: os.path.join(p, "bio/ngs/enrichment", "macs.rules")
if atac_cfg['trimadaptor']:
    include: os.path.join(p, "bio/ngs/qc", "cutadapt.rules")
if atac_cfg['bamfilter']:
    include: os.path.join(p, "bio/ngs/tools", "bamtools.rules")

    
ruleorder: picard_merge_sam > picard_sort_bam 
ruleorder: picard_sort_bam > picard_add_or_replace_read_groups
ruleorder: picard_add_or_replace_read_groups > picard_mark_duplicates
ruleorder: picard_mark_duplicates > atacseq_correct_coordinates
ruleorder: atacseq_correct_coordinates > bowtie_align
ruleorder: picard_sort_bam > bowtie_align
ruleorder: picard_merge_sam > bowtie_align
ruleorder: picard_mark_duplicates > bowtie_align

# Set temporary and protected outputs
set_output(workflow,
           temp_rules = main_cfg['temp_rules'] + main_cfg['temp_rules_default'],
           temp_filetypes=main_cfg['temp_filetypes'] + main_cfg['temp_filetypes_default'],
           protected_rules = main_cfg['protected_rules'] + main_cfg['protected_rules_default'],
           protected_filetypes=main_cfg['protected_filetypes'] + main_cfg['protected_filetypes_default'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

##############################
# Targets
##############################
ALIGN_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].run_id_re, 
    src_re = ngs_cfg['sampleorg'].raw_run_re, 
    target_suffix = ALIGN_TARGET_SUFFIX, 
    **ngs_cfg)

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].sample_re, 
    src_re = ngs_cfg['sampleorg'].raw_run_re, 
    target_suffix = MERGE_TARGET_SUFFIX, 
    **ngs_cfg)

PREFIX = ".sort.merge.filter" if atac_cfg['bamfilter'] else ".sort.merge"

DFILTER_TARGET_SUFFIX = PREFIX + ".offset.dfilt.bed"
DFILTER_TARGETS = []
if 'dfilter' in atac_cfg['peakcallers']:
    DFILTER_TARGETS = generic_target_generator(
        tgt_re = ngs_cfg['sampleorg'].sample_re, 
        src_re = ngs_cfg['sampleorg'].raw_run_re, 
        target_suffix = DFILTER_TARGET_SUFFIX, 
        **ngs_cfg) 

MACS2_TARGET_SUFFIX = PREFIX + ".offset_peaks.xls"
MACS2_TARGETS = []
if 'macs2' in atac_cfg['peakcallers']:
    MACS2_TARGETS = generic_target_generator(
        tgt_re = ngs_cfg['sampleorg'].sample_re, 
        src_re = ngs_cfg['sampleorg'].raw_run_re, 
        target_suffix = MACS2_TARGET_SUFFIX,  
        **ngs_cfg) 

DUP_METRICS_SUFFIX = ".sort.merge.dup.dup_metrics"
DUP_METRICS_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].sample_re, 
    target_suffix =  DUP_METRICS_SUFFIX, 
    src_re = ngs_cfg['sampleorg'].raw_run_re, 
    **ngs_cfg)

ALIGN_METRICS_SUFFIX = ".sort.merge.dup.align_metrics"
ALIGN_METRICS_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].sample_re, 
    target_suffix =  ALIGN_METRICS_SUFFIX, 
    src_re = ngs_cfg['sampleorg'].raw_run_re, 
    **ngs_cfg)

INSERT_METRICS_SUFFIX = ".sort.merge.dup.insert_metrics"
INSERT_METRICS_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].sample_re, 
    target_suffix =  INSERT_METRICS_SUFFIX, 
    src_re = ngs_cfg['sampleorg'].raw_run_re, 
    **ngs_cfg)

REPORT_TARGETS = ["report/atacseq_all_rulegraph.png", "report/atacseq_summary.html"]

BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in DFILTER_TARGETS] +\
                 [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in MACS2_TARGETS] +\
                 [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in MACS2_TARGETS]


# Rules
rule atacseq_all:
    """Run ATAC-seq pipeline"""
    input: DFILTER_TARGETS + MACS2_TARGETS + DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS + REPORT_TARGETS

rule atacseq_align:
    """Run ATAC-seq alignment"""
    input: ALIGN_TARGETS

rule atacseq_merge:
    """Run ATAC-seq alignment, duplication removal and merge"""
    input: MERGE_TARGETS

rule atacseq_metrics:
    """Run ATAC-seq alignment and corresponding metrics only"""
    input: DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS

rule atacseq_bigwig:
    """Convert peak-calling bed output to bigwig"""
    input: BIGWIG_TARGETS

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

rule atacseq_report:
    """Write report"""
    input: cutadapt = os.path.join("{path}", "cutadapt.summary.csv") if atac_cfg['trimadaptor'] else [],
           picard = [("report/picard.sort.merge.dup{sfx}.metrics.csv".format(sfx=sfx), 
           "report/picard.sort.merge.dup{sfx}.hist.csv".format(sfx=sfx)) for sfx in [workflow._rules[x].params.suffix for x in picard_config['qcrules']]],
           qualimap = [os.path.join("{path}", "sample{}.qualimap.globals.csv".format(MERGE_TARGET_SUFFIX)),
                       os.path.join("{path}", "sample{}.qualimap.coverage_per_contig.csv".format(MERGE_TARGET_SUFFIX))],
           rulegraph = "report/atacseq_all_rulegraph.png"
    output: html = os.path.join("{path}", "atacseq_summary.html")
    run:
        d = {}
        env = Environment(loader = PackageLoader("snakemakelib", "_templates"))
        tp = env.get_template('workflow_atacseq_qc.html')
        if atac_cfg['trimadaptor']:
            d.update({'cutadapt' : make_cutadapt_summary_plot(input.cutadapt)})
        d.update({'qualimap' : make_qualimap_plots(*input.qualimap)})
        d.update({'picard' : make_picard_summary_plots(input.picard)})
        d.update({'rulegraph' : {'uri' : data_uri(input.rulegraph), 'file' : input.rulegraph, 'fig' : input.rulegraph,
                                 'target' : 'atacseq_all'}})
        with open(output.html, "w") as fh:
            fh.write(static_html(tp, **d))

                
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

