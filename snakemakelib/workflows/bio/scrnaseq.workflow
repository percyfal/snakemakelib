# -*- snakemake -*-
import shutil
import os
from os.path import join, dirname, relpath, exists
import pickle
from snakemake.report import data_uri
from bokehutils.publish import static_html
from snakemake.workflow import workflow
from snakemakelib.io import set_output
from snakemakelib.utils import SmlTemplateEnv
from snakemakelib.config import update_config, SNAKEMAKELIB_RULES_PATH
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.workflow.scrnaseq import scrnaseq_alignment_qc_plots, scrnaseq_pca_plots, read_gene_expression, pca, pca_results, number_of_detected_genes

def _merge_suffix(aligner, quantification=[]):
    align_config = config['bio.ngs.align.' + aligner]
    if aligner == "star":
        return align_config['align']['suffix'].replace('.bam', '_unique.bam')
    elif aligner in ["bowtie", "bowtie2"]:
        return '_unique.bam'

def _merge_tx_suffix(aligner, quantification=[]):
    align_config = config['bio.ngs.align.' + aligner]
    if aligner == "star":
        return ".Aligned.toTranscriptome.out_unique.bam"

def _find_transcript_bam(wildcards):
    sources = generic_target_generator(
        tgt_re=config['bio.ngs.settings']['sampleorg'].run_id_re,
        target_suffix = _merge_tx_suffix(config['bio.ngs.settings']['aligner'], config['bio.ngs.settings']['rnaseq']['quantification']),
        **config['bio.ngs.settings'])
    sources = [src for src in sources if dirname(src).startswith(wildcards.prefix)]
    return sources

def find_scrnaseq_merge_inputs(wildcards):
    """Find platform unit specific aligned bam files as input to picard
    merge. 

    NB: these are *not* the transcript-specific alignment files.
    """
    sources = generic_target_generator(
        tgt_re=config['bio.ngs.settings']['sampleorg'].run_id_re, 
        target_suffix = _merge_suffix(config['bio.ngs.settings']['aligner'], config['bio.ngs.settings']['rnaseq']['quantification']), 
        **config['bio.ngs.settings'])
    sources = [src for src in sources if dirname(src).startswith(wildcards.prefix)]
    return sources


# Configuration
config_default = {
    'settings' : {
        'temp_rules' : [],
        'temp_rules_default' : ['sratools_prefetch', 'star_align', 'bamtools_filter_unique'],
    },
    'workflows.bio.scrnaseq' : {
        'qc' : {
            # Add parameters here
        },
        'quantification' :  ['rsem', 'rpkmforgenes'],
        'db' : {
            'do_multo' : False,  ## Set to True to run generate multo database
        },
        'metadata' : None, ## Sample metadata
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
    'bio.ngs.qc.rseqc' : {
        'rules' : ['rseqc_qc_8', 'rseqc_qc_8_summary'],
    },
    'bio.ngs.rnaseq.rpkmforgenes' : {
        'rules' : ['rpkmforgenes_from_bam', 'rpkmforgenes_summarize_expression_data'],
    },
    'bio.ngs.rnaseq.rsem' : {
        'rules' : ['rsem_calculate_expression', 'rsem_prepare_reference', 'rsem_summarize_expression_data'],
    },

}

update_config(config_default, config)
config = config_default

aligner = config['bio.ngs.settings']['aligner']

# Include necessary snakemakelib rules
include: join(SNAKEMAKELIB_RULES_PATH, 'settings.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'utils.rules')
align_include = join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/align", aligner + ".rules")
include: align_include
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/tools", "bamtools.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "rseqc.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "picard.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/tools", "samtools.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/rnaseq", "rpkmforgenes.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/rnaseq", "rsem.rules")

if aligner in ["bowtie", "bowtie2"]:
    ruleorder: bamtools_filter_unique > picard_merge_sam > picard_sort_bam > bowtie_align

if config['workflows.bio.scrnaseq']['db']['do_multo']:
    include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/tools", "multo.rules")
    
if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'scrnaseq.workflow'; set workdir before include statement!")

# Set temporary and protected outputs
set_output(workflow,
           temp_rules = config['settings']['temp_rules'] + config['settings']['temp_rules_default'],
           temp_filetypes=config['settings']['temp_filetypes'] + config['settings']['temp_filetypes_default'],
           protected_rules = config['settings']['protected_rules'] + config['settings']['protected_rules_default'],
           protected_filetypes=config['settings']['protected_filetypes'] + config['settings']['protected_filetypes_default'])

##################################################
# Target definitions
##################################################
ALIGN_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re,
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re,
    target_suffix = _merge_suffix(aligner),
    **config['bio.ngs.settings'])

RSEQC_TARGETS = generic_target_generator(
    tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix = '.merge_rseqc/rseqc_qc_8.txt',
    src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

RPKMFORGENES_TARGETS = []
if 'rpkmforgenes' in config['workflows.bio.scrnaseq']['quantification']:
    RPKMFORGENES_TARGETS = generic_target_generator(
        tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re,
        target_suffix = '.merge.rpkmforgenes',
        src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re,
        **config['bio.ngs.settings']) + ['report/rpkmforgenes.merge.csv']

RSEM_TARGETS = []
if 'rsem' in config['workflows.bio.scrnaseq']['quantification']:
    RSEM_TARGETS = generic_target_generator(
        tgt_re = config['bio.ngs.settings']['sampleorg'].sample_re,
        target_suffix = '.merge.tx.isoforms.results',
        src_re = config['bio.ngs.settings']['sampleorg'].raw_run_re,
        **config['bio.ngs.settings']) + ['report/rsem.merge.tx.genes.csv', 'report/rsem.merge.tx.isoforms.csv']

    
REPORT_TARGETS = ['report/star.Aligned.out.csv', 'report/scrnaseq_summary.html']

# Additional merge rule for transcript alignment files
rule scrnaseq_picard_merge_sam_transcript:
    """scrnaseq picard: merge sam files from transcript alignments.

    NB: always outputs bam files!
    """
    params: cmd = config['bio.ngs.qc.picard']['cmd'] + "MergeSamFiles",
            options = " ".join([config['bio.ngs.qc.picard']['options'],
                                config['bio.ngs.qc.picard']['merge_sam']['options']])
    input: _find_transcript_bam
    output: merge="{path}" + os.sep + "{prefix}." + "merge.tx.bam"
    run: 
      if (len(input) > 1):
          inputstr = " ".join(["INPUT={}".format(x) for x in input])
          shell("{cmd} {ips} OUTPUT={out} {opt}".format(cmd=params.cmd, ips=inputstr, out=output.merge, opt=params.options))
      else:
          if exists(output.merge):
              os.unlink(output.merge)
          shutil.copy(input[0], output.merge)

# QC rules
QC_INPUT = []
rule scrnaseq_qc:
    """Perform basic qc on samples"""
    input: starcsv = join("{path}", "star.Aligned.out.csv"),
           rseqc_read_distribution = join("{path}", "read_distribution_summary_merge_rseqc.csv"),
           rseqc_gene_coverage = join("{path}", "gene_coverage_summary_merge_rseqc.csv"),
           rsemgenes = join("{path}", "rsem.merge.tx.genes.csv") if 'rsem' in config['workflows.bio.scrnaseq']['quantification'] else [],
           rsemisoforms = join("{path}", "rsem.merge.tx.isoforms.csv")  if 'rsem' in config['workflows.bio.scrnaseq']['quantification'] else [],
           rsemgenespca = join("{path}", "rsem.merge.tx.genes.pca.csv") if 'rsem' in config['workflows.bio.scrnaseq']['quantification'] else [],
           rpkmforgenes = join("{path}", "rpkmforgenes.merge.csv") if 'rpkmforgenes' in config['workflows.bio.scrnaseq']['quantification'] else [],
           rpkmforgenespca = join("{path}", "rpkmforgenes.merge.pca.csv") if 'rpkmforgenes' in config['workflows.bio.scrnaseq']['quantification'] else [],
           rulegraph = join("{path}", "scrnaseq_all_rulegraph.png"),
           globalconf = join("{path}", "smlconf_global.yaml")
    output: html = join("{path}", "scrnaseq_summary.html")
    run:
        d = {'align': scrnaseq_alignment_qc_plots(input.rseqc_read_distribution,
                                                  input.rseqc_gene_coverage,
                                                  input.starcsv)}
        d['align'].update({
             'uri': [data_uri(input.rseqc_read_distribution),
                     data_uri(input.rseqc_gene_coverage),
                     data_uri(input.starcsv)],
             'file': [input.rseqc_read_distribution,
                     input.rseqc_gene_coverage,
                     input.starcsv],
        })
        d.update({'rulegraph' : {'uri' : data_uri(input.rulegraph), 'file' : input.rulegraph, 'fig' : input.rulegraph, 'target' : 'scrnaseq_all'}})
        if input.rsemgenes:
            #d.update({'rsem' : {'file': [input.rsemgenes, input.rsemisoforms]},
            d.update({'rsem' : {'file': [input.rsemgenespca], 'uri': [data_uri(input.rsemgenespca)]}})
            d['rsem'].update(scrnaseq_pca_plots(input.rsemgenespca,
                                                metadata=config['workflows.bio.scrnaseq']['metadata'],
                                                pcaobjfile=input.rsemgenespca.replace(".pca.csv",
                                                                                      ".pcaobj.pickle")))
        if input.rpkmforgenespca:
            d.update({'rpkmforgenes' : {'file': [input.rpkmforgenespca], 'uri': [data_uri(input.rpkmforgenespca)]}})
            d['rpkmforgenes'].update(
                scrnaseq_pca_plots(input.rpkmforgenespca,
                                   metadata=config['workflows.bio.scrnaseq']['metadata'],
                                   pcaobjfile=input.rpkmforgenespca.replace(".pca.csv", ".pcaobj.pickle")))
        d.update({'version' : config['_version'], 'config' : {'uri' : data_uri(input.globalconf), 'file' : input.globalconf}})
        tp = SmlTemplateEnv.get_template('workflow_scrnaseq_qc.html')
        with open(output.html, "w") as fh:
            fh.write(static_html(tp, **d))

rule scrnaseq_pca:
    input: expr = "{prefix}.csv",
           annotation = config['bio.ngs.settings']['annotation']['transcript_annot_gtf'] if config['bio.ngs.settings']['annotation']['transcript_annot_gtf'] else []
    output: pca = "{prefix}.pca.csv", pcaobj = "{prefix}.pcaobj.pickle"
    run:
        import math
        expr_long = read_gene_expression(input.expr,
                                         annotation=input.annotation)
        expr_long["TPM"] = [math.log2(x+1.0) for x in expr_long["TPM"]]
        detected_genes = number_of_detected_genes(expr_long)
        # FIXME: configurations?
        expr = expr_long.pivot_table(columns="gene_id", values="TPM",
                                     index="sample")
        pcaobj = pca(expr)
        pcares = pca_results(pcaobj, expr, metadata=config['workflows.bio.scrnaseq']['metadata'])
        if not detected_genes is None:
            pcares = pcares.join(detected_genes)
        with open(output.pca, "w") as fh:
            pcares.to_csv(fh)
        with open(output.pcaobj, "wb") as fh:
            pickle.dump(pcaobj, fh)
        
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
