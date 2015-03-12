# Copyright (C) 2014 by Per Unneberg
from snakemakelib.config import get_sml_config

sml_config = get_sml_config()

def report_label():
    """Return the report label based on bowtie2 and paired_end flags."""
    cfg = get_sml_config('bio.ngs.methylseq.bismark')
    report_label = "_bismark"
    report_label += "_bt2" if cfg['bowtie2'] else ""
    report_label += "_PE_report" if cfg['paired_end'] else "_SE_report"
    return report_label

def align_suffix(basename=False):
    """Return the align suffix based on bowtie2 and paired_end flags.
    
    Args:
      basename: return only the basename, stripping out file suffix and possible _pe label

    """
    cfg = get_sml_config('bio.ngs.methylseq.bismark')
    align_suffix = "_bismark"
    align_suffix += "_bt2" if cfg['bowtie2'] else ""
    if not basename:
        align_suffix += "_pe" if cfg['paired_end'] else ""
        align_suffix += ".bam"
    return align_suffix

