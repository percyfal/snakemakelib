# Copyright (C) 2014 by Per Unneberg

def report_label(bismark_cfg):
    """Return the report label based on bowtie2 and paired_end flags."""
    report_label = "_bismark"
    report_label += "_bt2" if bismark_cfg['bowtie2'] else ""
    report_label += "_PE_report" if bismark_cfg['paired_end'] else "_SE_report"
    return report_label

def align_suffix(bismark_cfg, basename=False):
    """Return the align suffix based on bowtie2 and paired_end flags.

    Args:
      bismark_cfg: configuration for 'bio.ngs.methylseq.bismark'
      basename: return only the basename, stripping out file suffix and possible _pe label

    """
    align_suffix = "_bismark"
    align_suffix += "_bt2" if bismark_cfg['bowtie2'] else ""
    if not basename:
        align_suffix += "_pe" if bismark_cfg['paired_end'] else ""
        align_suffix += ".bam"
    return align_suffix

