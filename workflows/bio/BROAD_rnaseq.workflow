# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config
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
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs", "settings.rules")
include: os.path.join(p, "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(p, "bio/ngs/qc", "picard.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "tuxedo.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "rnaseqc.rules")
include: os.path.join(p, "bio/ngs/rnaseq", "rsem.rules")

ruleorder: tuxedo_tophat > tuxedo_bowtie_align > samtools_sam2bam

# Get relevant config sections
cfg = get_sml_config('bio.ngs.settings')
tux_cfg = get_sml_config('bio.ngs.rnaseq.tuxedo')
path = cfg.get('path') if not cfg.get('path') is None else os.curdir

# Targets
RSEM_TARGETS = generic_target_generator(fmt=ngs_cfg['sample_pfx_fmt'] + '.merge.isoforms.results', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

TOPHAT2_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.tophat2', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

CUFFLINKS_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.cufflinks', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

CUFFLINKS_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.cufflinks_quant', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

RNASEQC_TARGETS = generic_target_generator(fmt=ngs_cfg['run_id_pfx_fmt'] + '.tophat2/accepted_hits.rg.resorted.dup.rnaseqc', rg=ReadGroup(ngs_cfg['run_id_pfx_re'] + ngs_cfg['read1_label'] + ngs_cfg['fastq_suffix']), cfg=ngs_cfg, path=path)

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
