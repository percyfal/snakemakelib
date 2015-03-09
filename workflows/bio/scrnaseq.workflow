# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, get_sml_config

scrnaseq_config = {
    
}

update_sml_config(scrnaseq_config)

# Include necessary snakemakelib rules
p = os.path.join(os.pardir, os.pardir, 'rules')
include: os.path.join(p, 'settings.rules')
include: os.path.join(p, 'utils.rules')
include: os.path.join(p, "bio/ngs/align", "star.rules")
include: os.path.join(p, "bio/ngs/qc", "rseqc.rules")
include: os.path.join(p, "bio/ngs/tools", "bamtools.rules")
include: os.path.join(p, "bio/ngs/tools", "samtools.rules")

# All rules
rule scrnaseq_all:
    """Run scRNAseq pipeline"""
    input: STAR_TARGETS + QC_TARGETS + RPKM_TARGETS
