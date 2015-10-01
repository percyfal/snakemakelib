# Copyright (C) 2015 by Per Unneberg
import math
import pandas as pd
from snakemakelib.log import LoggerManager

logger = LoggerManager().getLogger(__name__)

def number_of_detected_genes(expr, cutoff=1.0, **kwargs):
    """Aggregate expression data frame to count number of detected genes

    Args:
      expr (DataFrame): pandas data frame with expression values
      cutoff (float): cutoff for detected gene

    Returns:
      detected_genes (DataFrame): aggregated data fram with number of detected genes per sample
    """
    expr_long = read_gene_expression(expr)
    expr_long["TPM"] = [math.log2(x+1.0) for x in expr_long["TPM"]]
    try:
        detected_genes = expr_long.groupby("sample").agg(lambda x: sum(x > cutoff))
    except Exception as e:
        logger.warning("Failed to group genes by sample :", e)
        detected_genes = None
    return detected_genes

def _gene_name_map_from_gtf(gtf, unit_id, unit_name):
    """Get a mapping from gene_id to gene_name"""
    mapping = {}
    for feature in gtf[8]:
        tmp = {k.replace("\"", ""):v.replace("\"", "") for k, v in [x.split(" ") for x in feature.split("; ")]}
        mapping[tmp.get(unit_id, "")] = tmp.get(unit_name, tmp.get(unit_id, ""))
    return mapping


def read_gene_expression(infile, annotation=None, unit_id="gene_id",
                         unit_name="gene_name"):
    """Read gene expression file, renaming genes if annotation present.

    NB: currently assumes annotation file is in gtf format and that
    gene expression levels, not transcript, are used

    Args:
      infile (str): infile name
      annotation (str): annotation file, gtf format
      unit_id (str): id of measurement unit; gene_id or transcript_id
      unit_name (str): name of measurement unit, as defined by annotation file
    
    Returns:
      expr (DataFrame): (possibly annotated) data frame

    """
    expr = pd.read_csv(infile)
    if annotation:
        annot = pd.read_table(annotation, header=None)
        mapping = _gene_name_map_from_gtf(annot, unit_id, unit_name)
        expr[unit_name] = expr[unit_id].map(mapping.get)
    return expr
