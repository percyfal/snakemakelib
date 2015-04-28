# Copyright (C) 2015 by Per Unneberg
import os
from snakemake.workflow import config
from snakemakelib.config import get_sml_config
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


def ref():
    """Return the fasta reference sequence as a string
    
    Args:
      None

    Returns:
      seq: A <string> representing the reference fasta sequence
    """
    cfg = config["bio.ngs.settings"]
    # If no build issue a warning
    if not cfg['db']['build']:
        raise ValueError('No build defined: cannot generate reference fasta file name. Either provide a build, a reference file name or write a custom ref() function that returns the full path name of the reference fasta file.')
    # Otherwise, see if build_config is defined
    if not cfg['db']['build_config']:
        raise ValueError('build defined but no build_config: cannot generate reference fasta file name. Either provide a build_config file in which the build defines reference and index files or write a custom ref() function that returns the full path name to the reference fasta file.')
    # TODO: make sure build_config points to a file; read file and parse to find reference
    return None
    
def index(application):
    """Return the index files for a given application.

    Args:
      application: a <string> representing the application.

    Returns:
      index: a <string> representing the index in the format required by the particular application
    """
    ngs_cfg = config["bio.ngs.settings"]

    if not ngs_cfg['db']['build_config']:
        smllogger.debug("No build_config present: assuming index locations are organized according to cloudbiolinux conventions")
        if application in ["bwa", "rsem", "bowtie"]:
            prefix, _ = os.path.splitext(ngs_cfg['db']['ref'])
        else:
            prefix = ngs_cfg['db']['ref']
        index_root = os.path.dirname(os.path.dirname(prefix))
        basename = os.path.basename(prefix)
        if application in ["star"]:
            return os.path.join(index_root, application, config['bio.ngs.align.star']['star_index']['genome'])
        else:
            return os.path.join(index_root, application, basename)
    else:
        return None
