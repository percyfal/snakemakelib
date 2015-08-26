# -*- snakemake -*-
import os
from bokehutils.publish import static_html
from snakemake.workflow import workflow
from snakemakelib.io import set_output
from snakemakelib.utils import SmlTemplateEnv
from snakemakelib.config import update_config
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.workflow.scrnaseq import scrnaseq_qc_plots

def _merge_suffix(aligner, quantification=[]):
    align_cfg = config['bio.ngs.align.' + aligner]
    if aligner == "star":
        return align_cfg['align']['suffix'].replace('.bam', '_unique.bam')
    elif aligner in ["bowtie", "bowtie2"]:
        return '_unique.bam'

def _merge_tx_suffix(aligner, quantification=[]):
    align_cfg = config['bio.ngs.align.' + aligner]
    if aligner == "star":
        return ".Aligned.toTranscriptome.out_unique.bam"

def _find_transcript_bam(wildcards):
    ngs_cfg = config['bio.ngs.settings']
    picard_cfg = config['bio.ngs.qc.picard']
    sources = generic_target_generator(
        tgt_re=ngs_cfg['sampleorg'].run_id_re,
        target_suffix = _merge_tx_suffix(ngs_cfg['aligner'], ngs_cfg['rnaseq']['quantification']),
        **ngs_cfg)
    sources = [src for src in sources if os.path.dirname(src).startswith(wildcards.prefix)]
    return sources

def find_scrnaseq_merge_inputs(wildcards):
    """Find platform unit specific aligned bam files as input to picard
    merge. 

    NB: these are *not* the transcript-specific alignment files.
    """
    ngs_cfg = config['bio.ngs.settings']
    sources = generic_target_generator(
        tgt_re=ngs_cfg['sampleorg'].run_id_re, 
        target_suffix = _merge_suffix(ngs_cfg['aligner'], ngs_cfg['rnaseq']['quantification']), 
        **ngs_cfg)
    sources = [src for src in sources if os.path.dirname(src).startswith(wildcards.prefix)]
    return sources


# Configuration
config_default = {
    'settings' : {
        'temp_rules' : [],
        'temp_rules_default' : ['sratools_prefetch', 'star_align', 'bamtools_filter'],
    },
    'workflows.bio.scrnaseq' : {
        'qc' : {
            
        },
        'quantification' :  ['rsem', 'rpkmforgenes']
    },
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

config = update_config(config_default, config)
main_cfg = config['settings']
ngs_cfg = config['bio.ngs.settings']
aligner = ngs_cfg['aligner']
work_cfg = config['workflows.bio.scrnaseq']

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

# Set temporary and protected outputs
set_output(workflow,
           temp_rules = main_cfg['temp_rules'] + main_cfg['temp_rules_default'],
           temp_filetypes=main_cfg['temp_filetypes'] + main_cfg['temp_filetypes_default'],
           protected_rules = main_cfg['protected_rules'] + main_cfg['protected_rules_default'],
           protected_filetypes=main_cfg['protected_filetypes'] + main_cfg['protected_filetypes_default'])

##################################################
# Target definitions
##################################################
ALIGN_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].run_id_re,
    src_re = ngs_cfg['sampleorg'].raw_run_re,
    target_suffix = _merge_suffix(aligner),
    **ngs_cfg)

RSEQC_TARGETS = generic_target_generator(
    tgt_re = ngs_cfg['sampleorg'].sample_re,
    target_suffix = '.merge_rseqc/rseqc_qc_8.txt',
    src_re = ngs_cfg['sampleorg'].raw_run_re,
    **ngs_cfg)

RPKMFORGENES_TARGETS = []
if 'rpkmforgenes' in work_cfg['quantification']:
    RPKMFORGENES_TARGETS = generic_target_generator(
        tgt_re = ngs_cfg['sampleorg'].sample_re,
        target_suffix = '.merge.rpkmforgenes',
        src_re = ngs_cfg['sampleorg'].raw_run_re,
        **ngs_cfg) 

RSEM_TARGETS = []
if 'rsem' in work_cfg['quantification']:
    RSEM_TARGETS = generic_target_generator(
        tgt_re = ngs_cfg['sampleorg'].sample_re,
        target_suffix = '.merge.tx.isoforms.results',
        src_re = ngs_cfg['sampleorg'].raw_run_re,
        **ngs_cfg) + ['report/rsem.merge.tx.genes.csv', 'report/rsem.merge.tx.isoforms.csv']

REPORT_TARGETS = ['report/star.Aligned.out.csv', 'report/star.Aligned.out.mapping_summary.html']

picard_config = config["bio.ngs.qc.picard"]

# Additional merge rule for transcript alignment files
rule scrnaseq_picard_merge_sam_transcript:
    """scrnaseq picard: merge sam files from transcript alignments.

    NB: always outputs bam files!
    """
    params: cmd = picard_config['cmd'] + "MergeSamFiles",
            options = " ".join([picard_config['options'],
                                picard_config['merge_sam']['options']])
    input: _find_transcript_bam
    output: merge="{path}" + os.sep + "{prefix}." + "merge.tx.bam"
    run: 
      if (len(input) > 1):
          inputstr = " ".join(["INPUT={}".format(x) for x in input])
          shell("{cmd} {ips} OUTPUT={out} {opt}".format(cmd=params.cmd, ips=inputstr, out=output.merge, opt=params.options))
      else:
          os.symlink(os.path.relpath(input[0], wildcards.path), output.merge)

# QC rules
QC_INPUT = []
rule scrnaseq_qc:
    """Perform basic qc on samples"""
    input: starcsv = os.path.join("{path}", "star.Aligned.out.csv"),
           rseqc_read_distribution = os.path.join("{path}", "read_distribution_summary_merge_rseqc.csv"),
           rseqc_gene_coverage = os.path.join("{path}", "gene_coverage_summary_merge_rseqc.csv"),
           rsemgenes = os.path.join("{path}", "rsem.merge.tx.genes.csv"),
           rsemisoforms = os.path.join("{path}", "rsem.merge.tx.isoforms.csv"),
           rulegraph = os.path.join("{path}", "scrnaseq_all.png"),
           globalconf = os.path.join("{path}", "smlconf_global.yaml")
    output: html = os.path.join("{path}", "scrnaseq_summary.html")
    run:
        d = {'align': scrnaseq_qc_plots(input.rseqc_read_distribution,
                                        input.rseqc_gene_coverage,
                                        input.starcsv)}
        d.update({'rulegraph' : {'uri' : data_uri(input.rulegraph), 'file' : input.rulegraph, 'fig' : input.rulegraph, 'target' : 'scrnaseq_all'}})
        d.update({'rsem' : {'file': [input.rsemgenes, input.rsemisoforms]}})

        d.update({'version' : config['_version'], 'config' : {'uri' : data_uri(input.globalconf), 'file' : input.globalconf}})
        tp = SmlTemplateEnv.get_template('workflow_scrnaseq_qc.html')
        with open(output.html, "w") as fh:
            fh.write(static_html(tp, **d))

rule scrnaseq_pca:
    input: csv = os.path.join("{path}", "rsem.merge.tx.{type}.csv")
    output: tmp = os.path.join("{path}", "rsem.merge.tx.{type}.pca")
    run:
        import pandas as pd
        from sklearn.decomposition import PCA
        #df = pd.read_csv(input.csv, index_col=0)
        #pca = PCA(n_components='mle')
        #pca.fit(df.head())
        #print (df.head())

# All rules
rule scrnaseq_all:
    """Run scRNAseq pipeline"""
    input: ALIGN_TARGETS + RSEQC_TARGETS + RPKMFORGENES_TARGETS + RSEM_TARGETS + REPORT_TARGETS

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

rule scrnaseq_clean:
    """Clean working directory. WARNING: will remove all files except
    (.fastq|.fastq.gz) and csv files
    """
    params: d = workflow._workdir
    shell: 'for f in `find  {params.d} -type f -name "*" | grep -v ".fastq$" | grep -v ".fastq.gz$" | grep -v ".csv$"`; do echo removing $f; rm -f $f; done'
