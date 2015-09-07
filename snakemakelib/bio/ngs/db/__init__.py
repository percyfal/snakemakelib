# Copyright (C) 2015 by Per Unneberg
import os
from snakemakelib.log import LoggerManager

def index(ref, application, build=None, index="", index_name=""):
    return index

def annotation(db_config, annotation='ref-transcripts.gtf', ignore_extra_ref=False, fmt="gtf"):
    return annotation

def ref(ref, db_config):
    return ref

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

