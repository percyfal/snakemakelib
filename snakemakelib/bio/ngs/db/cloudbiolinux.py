# Copyright (C) 2015 by Per Unneberg
import os
import snakemakelib.bio.ngs.db
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def annotation(db_config, annotation="ref-transcripts.gtf", ignore_extra_ref=False, fmt="gtf"):
    """Return the annotation as a string"""
    if os.path.isabs(annotation):
        if not annotation.endswith(fmt):
            (root, ext) = os.path.splitext(annotation)
            annotation = "{root}.{ext}".format(root=root, ext=fmt)
        return annotation
    if db_config['ref']:
        smllogger.debug("reference set: assuming index locations are organized according to cloudbiolinux conventions")
        if db_config['extra_ref'] and not ignore_extra_ref:
            smllogger.debug("extra references set: renaming annotation file")
            annotation = "ref-transcripts-" + "-".join([os.path.splitext(os.path.basename(x))[0] for x in db_config['extra_ref']]) + "." + fmt
        annotation = os.path.join(os.path.dirname(db_config['ref']), os.pardir, 'rnaseq', annotation)
        return annotation

def ref(ref, db_config):
    """Return the fasta reference sequence as a string
    
    Args:
      db_config: configuration object for section bio.ngs.settings.db

    Returns:
      seq: A <string> representing the reference fasta sequence
    """
    if os.path.isabs(ref):
        return ref
    if db_config['ref']:
        smllogger.debug("reference set: assuming index locations are organized according to cloudbiolinux conventions")
        ref = os.path.join(os.path.dirname(db_config['ref']), ref)
        return ref
    else:    
        # If no build issue a warning
        if not db_config['build']:
            raise ValueError('No build defined: cannot generate reference fasta file name. Either provide a build, a reference file name or write a custom ref() function that returns the full path name of the reference fasta file.')
            # Otherwise, see if build_config is defined
        elif not db_config['build_config']:
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

## Set the functions
smllogger.info("Setting snakemakelib.bio.ngs.db functions index, ref, and annotation to to cloudbiolinux variants")
snakemakelib.bio.ngs.db.index = index
snakemakelib.bio.ngs.db.ref = ref
snakemakelib.bio.ngs.db.annotation = annotation
