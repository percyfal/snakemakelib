# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import gzip
import io
import unittest
import logging
from subprocess import Popen, PIPE
from nose.tools import raises

logger = logging.getLogger(__name__)

hg19_dirs = ['bowtie', 'bwa', 'seq']

def install_ucsc_genome(build, chr, start=None, end=None):
    logger.info("Installing genome build {}, chr {}".format(build, chr))
    if start and end:
        logger.info("Installing region {}-{}".format(start, end))
    # Assume 50 characters per row
    lines = abs(int((end - start) / 50  )) + 1
    ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/chromosomes/".format(build)
    url = os.path.join(ucsc, "{}.fa.gz".format(chr))
    genomedir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'seq')
    if not os.path.exists(genomedir):
        logger.info("Creating {}".format(genomedir))
        os.makedirs(genomedir)
    try:
        outfile = os.path.join(genomedir, os.path.basename(url))
        if not os.path.exists(outfile.replace(".gz", "")):
            logger.info("Downloading {} from {} with curl".format(outfile, url))
            cl = ["curl", url, "-o", outfile + ".tmp"]
            p = Popen(cl, stdout=PIPE)
            stdout, stderr = p.communicate()
            if p.returncode == 0:
               logger.info("{} generation completed; renaming to {}".format(outfile + ".tmp", outfile))
               os.rename(outfile + ".tmp", outfile)
            else:
                raise Exception("{} failed: \n{}".format(cmd, " ".join([stderr])))
            cmd = 'zcat {} | head -{} > {}'.format(outfile, lines, outfile.replace(".gz", ""))
            pipe = Popen(['/bin/bash', '-c', cmd])
            stdout, stderr = pipe.communicate()
            if pipe.returncode == 0:
                logger.info("{} generation completed; removing {}".format(outfile.replace(".gz", ""), outfile))
                os.unlink(outfile)
            else:
                raise Exception("{} failed: \n{}".format(cmd, " ".join([stderr])))
    except:
        logger.info("Failed to download {}".format(os.path.join(genomedir, os.path.basename(url))))
    # Add links to bwa, bowtie, bismark
    for d in ['bwa', 'bowtie', 'bismark']:
        if not os.path.exists(os.path.join(genomedir, os.pardir, d)):
            os.mkdir(os.path.join(genomedir, os.pardir, d))
        if not os.path.exists(os.path.join(genomedir, os.pardir, d, os.path.basename(outfile).replace(".gz", ""))):
            os.symlink(outfile.replace(".gz", ""), os.path.join(genomedir, os.pardir, d, os.path.basename(outfile).replace(".gz", "")))
    return

def setUp():
    """Set up test fixtures"""
    logger.info("Setting up test fixtures")
    install_ucsc_genome(build="hg19", chr="chr11", start=0, end=2000000)
    
