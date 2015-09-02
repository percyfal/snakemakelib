# -*- snakemake -*-
from os.path import join
from snakemakelib.config import update_config, SNAKEMAKELIB_RULES_PATH
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.utils import ReadGroup

# Default configuration settings
rnaseq_config = {
    'bio.ngs.settings' : {
        'fastq_suffix' : ".P.qtrim.fq",
    },
}

# Start by including the general snakefile
# Include necessary snakemakelib rules
include: join(SNAKEMAKELIB_RULES_PATH, 'settings.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs", "settings.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "sequenceprocessing.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/qc", "picard.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/rnaseq", "tuxedo.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/rnaseq", "rnaseqc.rules")
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/rnaseq", "rsem.rules")

ruleorder: tuxedo_tophat > tuxedo_bowtie_align > samtools_sam2bam

# Targets
RSEM_TARGETS = generic_target_generator(fmt=config['bio.ngs.settings']['sample_pfx_fmt'] + '.merge.isoforms.results', rg=ReadGroup(config['bio.ngs.settings']['run_id_pfx_re'] + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix']), cfg=config['bio.ngs.settings'])

TOPHAT2_TARGETS = generic_target_generator(fmt=config['bio.ngs.settings']['run_id_pfx_fmt'] + '.tophat2', rg=ReadGroup(config['bio.ngs.settings']['run_id_pfx_re'] + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix']), cfg=config['bio.ngs.settings'])

CUFFLINKS_TARGETS = generic_target_generator(fmt=config['bio.ngs.settings']['run_id_pfx_fmt'] + '.cufflinks', rg=ReadGroup(config['bio.ngs.settings']['run_id_pfx_re'] + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix']), cfg=config['bio.ngs.settings'])

CUFFLINKS_TARGETS = generic_target_generator(fmt=config['bio.ngs.settings']['run_id_pfx_fmt'] + '.cufflinks_quant', rg=ReadGroup(config['bio.ngs.settings']['run_id_pfx_re'] + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix']), cfg=config['bio.ngs.settings'])

RNASEQC_TARGETS = generic_target_generator(fmt=config['bio.ngs.settings']['run_id_pfx_fmt'] + '.tophat2/accepted_hits.rg.resorted.dup.rnaseqc', rg=ReadGroup(config['bio.ngs.settings']['run_id_pfx_re'] + config['bio.ngs.settings']['read1_label'] + config['bio.ngs.settings']['fastq_suffix']), cfg=config['bio.ngs.settings'])

rule BROAD_rnaseq_all:
    """Run all the analyses"""
    input: RSEM_TARGETS + TOPHAT2_TARGETS + CUFFLINKS_TARGETS + CUFFLINKS_QUANT_TARGETS + RNASEQC_TARGETS

rule BROAD_rnaseq_rsem:
    input: RSEM_TARGETS

rule BROAD_rnaseq_tophat2:
    input: TOPHAT2_TARGETS

rule BROAD_rnaseq_cufflinks:
    input: CUFFLINKS_TARGETS

rule BROAD_rnaseq_cufflinks_quant:
    input: CUFFLINKS_QUANT_TARGETS

rule BROAD_rnaseq_rnaseqc:
    input: RNASEQC_TARGETS

# # rule rule_12:
# # 	input: " $(TOPHAT2_TARGETS)"
# # 	output: "{prefix}/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.bam"
# # 	shell: "echo "Checkpoint that all tophat2 targets have been completed before running rnaseqc""

rule list_targets:
    """List currently defined targets"""
    run:
      print ("List targets")
      print ("RSEM targets: ", RSEM_TARGETS)
      print ("RNAseQC targets: ", RNASEQC_TARGETS)
      print ("Cufflinks targets: ", CUFFLINKS_TARGETS)
