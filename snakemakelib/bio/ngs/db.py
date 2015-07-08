# Copyright (C) 2015 by Per Unneberg
import os
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def chromosomes(ref):
    """Return the chromosome names for a given reference.

    Currently assumes the existence of a dictionary file.
    """
    dictfile = ref.replace(".fa", ".dict")
    try:
        with open(dictfile, "r") as fh:
            chroms = [x.split()[1].replace("SN:", "") for x in fh.readlines() if x.startswith("@SQ")]
    except:
        raise Exception("failed to read dict file {dict}".format(dictfile))
    return chroms


def annotation(db_cfg, annotation="ref-transcripts.gtf", ignore_extra_ref=False, fmt="gtf"):
    """Return the annotation as a string"""
    if os.path.isabs(annotation):
        if not annotation.endswith(fmt):
            (root, ext) = os.path.splitext(annotation)
            annotation = "{root}.{ext}".format(root=root, ext=fmt)
        return annotation
    if db_cfg['ref']:
        smllogger.debug("reference set: assuming index locations are organized according to cloudbiolinux conventions")
        if db_cfg['extra_ref'] and not ignore_extra_ref:
            smllogger.debug("extra references set: renaming annotation file")
            annotation = "ref-transcripts-" + "-".join([os.path.splitext(os.path.basename(x))[0] for x in db_cfg['extra_ref']]) + "." + fmt
        annotation = os.path.join(os.path.dirname(db_cfg['ref']), os.pardir, 'rnaseq', annotation)
        return annotation

def ref(ref, db_cfg):
    """Return the fasta reference sequence as a string
    
    Args:
      db_cfg: configuration object for section bio.ngs.settings.db

    Returns:
      seq: A <string> representing the reference fasta sequence
    """
    if os.path.isabs(ref):
        return ref
    if db_cfg['ref']:
        smllogger.debug("reference set: assuming index locations are organized according to cloudbiolinux conventions")
        ref = os.path.join(os.path.dirname(db_cfg['ref']), ref)
        return ref
    else:    
        # If no build issue a warning
        if not db_cfg['build']:
            raise ValueError('No build defined: cannot generate reference fasta file name. Either provide a build, a reference file name or write a custom ref() function that returns the full path name of the reference fasta file.')
            # Otherwise, see if build_config is defined
        elif not db_cfg['build_config']:
            raise ValueError('build defined but no build_config: cannot generate reference fasta file name. Either provide a build_config file in which the build defines reference and index files or write a custom ref() function that returns the full path name to the reference fasta file.')
            # TODO: make sure build_config points to a file; read file and parse to find reference
    return None

def index(ref, application, build=None, index="", index_name=""):
    if os.path.isabs(index):
        return index
    else:
        smllogger.debug("index not an absolute path: assuming index locations are organized according to cloudbiolinux conventions")
    if not build:
        smllogger.debug("No build_config present: assuming index locations are organized according to cloudbiolinux conventions")
        if application in ["rsem", "bowtie", "bowtie2"]:
            prefix, _ = os.path.splitext(ref)
        else:
            prefix = ref
        if application in ["rsem", "star"]:
            subdir = "rnaseq"
        else:
            subdir = ""
        index_root = os.path.join(os.path.dirname(prefix), os.pardir, subdir, application, os.path.dirname(index))
        index_name = os.path.basename(prefix) if not index_name else index_name
        return os.path.join(index_root, index_name)
    else:
        return None

