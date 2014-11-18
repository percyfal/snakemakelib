# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, sml_rules_path, get_sml_config

methylation_config = {
        'bio.ngs.methylseq.bismark' : {
        'methXtract' : {
            'options' : "--ignore_r2 2 --counts  --gzip -p --no_overlap",
        },
    },
}

update_sml_config(methylation_config)

include: os.path.join(sml_rules_path(), 'settings.rules')
include: os.path.join(sml_rules_path(), 'utils.rules')
include: os.path.join(sml_rules_path(), "bio/ngs", "settings.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/methylseq", "bismark.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/qc", "sequenceprocessing.rules")
include: os.path.join(sml_rules_path(), "bio/ngs/qc", "picard.rules")

bismark_cfg = get_sml_config('bio.ngs.methylseq.bismark')
cfg = get_sml_config('bio.ngs.settings')

path = cfg.get('path') if not cfg.get('path') is None else os.curdir

FASTQC_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_1_fastqc.html {path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_2_fastqc.html".split(), sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

#BISMARK_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_bismark_bt2_pe.deduplicated.cov", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

BISMARK_TARGETS = expand("{path}/{sample}/{sample}_bismark_bt2_pe.merge.deduplicated.cov", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

#BISMARK_REPORT_TARGETS = expand("{path}/{sample}/{flowcell}/{lane}_{flowcell}_{sample}_bismark_bt2_PE_report.html", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

BISMARK_REPORT_TARGETS = expand("{path}/{sample}/{sample}_bismark_bt2_PE_report.html", sample=cfg['samples'], flowcell=cfg['flowcells'], lane=cfg['lanes'], path=path)

# All rules
rule all:
    """Run all the analyses"""
    input: FASTQC_TARGETS + BISMARK_TARGETS + BISMARK_REPORT_TARGETS

rule run_fastqc:
    """Fastqc target rule. Run fastqc on files defined in FASTQC_TARGETS"""
    input: FASTQC_TARGETS

rule run_bismark:
    """bismark target rule. Run bismark on files defined in BISMARK_TARGETS"""
    input: BISMARK_TARGETS

# TODO: generic rule
rule list_targets:
    """List currently defined targets"""
    run:
      print (SAMPLE_TARGETS)
      print (FASTQC_TARGETS)
      print (BISMARK_TARGETS)

