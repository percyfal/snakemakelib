# -*- snakemake -*-
from os.path import join
import pysam
from jinja2 import Environment, PackageLoader
from bokehutils.publish import static_html
from snakemake.report import data_uri
from snakemake.utils import update_config
from snakemakelib.io import set_output
from snakemakelib.config import SNAKEMAKELIB_RULES_PATH
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.qc.cutadapt import make_cutadapt_summary_plot
from snakemakelib.bio.ngs.qc.qualimap import make_qualimap_plots
from snakemakelib.bio.ngs.qc.picard import make_picard_summary_plots

##############################
# Functions
##############################
def _merge_suffix():
    """Determine the merge suffix of the run files"""
    if config['workflows.bio.atac_seq']['trimadaptor']:
        return ".trimmed.sort.bam"
    else:
        return ".sort.bam"


##############################
# Default configuration settings custom-tailored for ATAC-Seq analysis
##############################
atac_config = {
    'workflows.bio.atac_seq' : {
        'aligner' : 'bowtie2',
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
    'bio.ngs.qc.cutadapt' : {
        'rules' : ['cutadapt_cut_paired_end', 'cutadapt_qc_summary'],
    },
    'bio.ngs.db.ucsc' : {
        'rules' : ['ucsc_wigToBigWig', 'ucsc_bedGraphToBigWig', 'ucsc_fetchChromSizes'],
    },
}
aligner_config = {
    'bio.ngs.align.bowtie' : {
        'align' : {
            'options' : '-X 2000',
        },
    },
    'bio.ngs.align.bowtie2' : {
        'align' : {
            'options' : '-X 2000',
        },
    },
    'bio.ngs.align.bwa' : {
    },
}

update_config(atac_config, config)
config = atac_config

ALIGN_TARGET_SUFFIX = ".bam"
if config['workflows.bio.atac_seq']['trimadaptor']:
    config['bio.ngs.qc.picard']['merge_sam']['suffix'] = '.trimmed.sort.bam'
    ALIGN_TARGET_SUFFIX = ".trimmed.bam"

aligner = config['workflows.bio.atac_seq']['aligner']
key = 'bio.ngs.align.' + aligner

config_default = aligner_config[key]
update_config(config_default, config)
config = config_default

##############################
# Include statements
##############################
include: join(SNAKEMAKELIB_RULES_PATH, 'settings.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs", "settings.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/align", aligner + ".rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/align", "blat.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/db", "ucsc.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "picard.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "qualimap.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/chromatin", "danpos.rules")
if 'dfilter' in config['workflows.bio.atac_seq']['peakcallers']:
    include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/enrichment", "dfilter.rules")
if 'macs2' in config['workflows.bio.atac_seq']['peakcallers']:
    include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/enrichment", "macs2.rules")
if config['workflows.bio.atac_seq']['trimadaptor']:
    include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "cutadapt.rules")
if config['workflows.bio.atac_seq']['bamfilter']:
    include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/tools", "bamtools.rules")

    
ruleorder: picard_merge_sam > picard_sort_sam 
ruleorder: picard_sort_sam > picard_add_or_replace_read_groups
ruleorder: picard_add_or_replace_read_groups > picard_mark_duplicates
ruleorder: picard_mark_duplicates > atacseq_correct_coordinates

# Set temporary and protected outputs
set_output(workflow,
           temp_rules = config['settings']['temp_rules'] + config['settings']['temp_rules_default'],
           temp_filetypes=config['settings']['temp_filetypes'] + config['settings']['temp_filetypes_default'],
           protected_rules = config['settings']['protected_rules'] + config['settings']['protected_rules_default'],
           protected_filetypes=config['settings']['protected_filetypes'] + config['settings']['protected_filetypes_default'])

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'ATAC-seq.workflow'; set workdir before include statement!")

##############################
# Targets
##############################
ALIGN_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
    target_suffix = ALIGN_TARGET_SUFFIX, 
    **config['bio.ngs.settings'])

MERGE_TARGET_SUFFIX = ".sort.merge.bam"
MERGE_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
    target_suffix = MERGE_TARGET_SUFFIX, 
    **config['bio.ngs.settings'])

PREFIX = ".sort.merge.filter" if config['workflows.bio.atac_seq']['bamfilter'] else ".sort.merge"

DFILTER_TARGET_SUFFIX = PREFIX + ".offset.dfilt.bed"
DFILTER_TARGETS = []
if 'dfilter' in config['workflows.bio.atac_seq']['peakcallers']:
    DFILTER_TARGETS = generic_target_generator(
        tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
        src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
        target_suffix = DFILTER_TARGET_SUFFIX, 
        **config['bio.ngs.settings']) 

MACS2_TARGET_SUFFIX = PREFIX + ".offset_peaks.xls"
MACS2_TARGETS = []
if 'macs2' in config['workflows.bio.atac_seq']['peakcallers']:
    MACS2_TARGETS = generic_target_generator(
        tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
        src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
        target_suffix = MACS2_TARGET_SUFFIX,  
        **config['bio.ngs.settings']) 

DUP_METRICS_SUFFIX = ".sort.merge.dup.dup_metrics"
DUP_METRICS_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
    target_suffix =  DUP_METRICS_SUFFIX, 
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
    **config['bio.ngs.settings'])

ALIGN_METRICS_SUFFIX = ".sort.merge.dup.align_metrics"
ALIGN_METRICS_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
    target_suffix =  ALIGN_METRICS_SUFFIX, 
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
    **config['bio.ngs.settings'])

INSERT_METRICS_SUFFIX = ".sort.merge.dup.insert_metrics"
INSERT_METRICS_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re, 
    target_suffix =  INSERT_METRICS_SUFFIX, 
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re, 
    **config['bio.ngs.settings'])

REPORT_TARGETS = ["report/atacseq_all_rulegraph.png", "report/atacseq_summary.html"]

BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in DFILTER_TARGETS] +\
                 [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in MACS2_TARGETS] +\
                 [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in MACS2_TARGETS]

##############################
# Collection rules
##############################
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


##############################
# Specific rules
##############################    
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
    input: cutadapt = join("{path}", "cutadapt.summary.csv") if config['workflows.bio.atac_seq']['trimadaptor'] else [],
           picard = [("report/picard.sort.merge.dup{sfx}.metrics.csv".format(sfx=sfx), 
           "report/picard.sort.merge.dup{sfx}.hist.csv".format(sfx=sfx)) for sfx in [workflow._rules[x].params.suffix for x in config['bio.ngs.qc.picard']['qcrules']]],
           qualimap = [join("{path}", "sample{}.qualimap.globals.csv".format(MERGE_TARGET_SUFFIX)),
                       join("{path}", "sample{}.qualimap.coverage_per_contig.csv".format(MERGE_TARGET_SUFFIX))],
           rulegraph = "report/atacseq_all_rulegraph.png"
    output: html = join("{path}", "atacseq_summary.html")
    run:
        d = {}
        env = Environment(loader = PackageLoader("snakemakelib", "_templates"))
        tp = env.get_template('workflow_atacseq_qc.html')
        if config['workflows.bio.atac_seq']['trimadaptor']:
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

