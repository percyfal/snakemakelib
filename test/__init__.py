# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import gzip
import io
import unittest
import logging
import time
from subprocess import Popen, PIPE
from nose.tools import raises

logger = logging.getLogger(__name__)

def install_align_ref():
    genomedir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', "Hsapiens", "hg19", 'seq')
    chr = "chr11"
    outfile = os.path.join(genomedir, "{}.fa.gz".format(chr))
    # Add links to bwa, bowtie
    for d in ['bwa', 'bowtie']:
        p = os.path.join(genomedir, os.pardir, d)
        if not os.path.exists(p):
            logger.info("Creating application directory {}".format(p))
            os.mkdir(p)
        p = os.path.join(genomedir, os.pardir, d, os.path.basename(outfile).replace(".gz", ""))
        if not os.path.exists(p):
            logger.info("Creating application file {}".format(p))
            os.symlink(outfile.replace(".gz", ""), p)
    return


def setUp():
    """Set up test fixtures"""
    logger.info("Setting up test fixtures")
    install_align_ref()
    
