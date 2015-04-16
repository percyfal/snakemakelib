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

config = {
    'hg19' : {
        'chr' : 'chr11',
        'start' : 0,
        'end' : 2000000,
        'species' : 'Hsapiens',

        'annotation' : {'url': 'ftp://ftp.sanger.ac.uk/pub/gencode/release_7/gencode.v7.annotation.gtf.gz',
                        'dest' : 'gencode.v7.annotation.chr11.gtf',
                        'lines' : 10000,
                    },
    },
}


def install_ucsc_genome(build, start=None, end=None):
    chr = config[build]['chr']
    logger.info("Installing genome build {}, chr {}".format(build, chr))
    if start and end:
        logger.info("Installing region {}-{}".format(start, end))
    # Assume 50 characters per row
    lines = abs(int((end - start) / 50  )) + 1
    ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/chromosomes/".format(build)
    url = os.path.join(ucsc, "{}.fa.gz".format(chr))
    genomedir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', config[build]['species'], build, 'seq')
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
    for d in ['bwa', 'bowtie', 'bismark', 'star']:
        p = os.path.join(genomedir, os.pardir, d)
        if not os.path.exists(p):
            logger.info("Creating application directory {}".format(p))
            os.mkdir(p)
        p = os.path.join(genomedir, os.pardir, d, os.path.basename(outfile).replace(".gz", ""))
        if not os.path.exists(p):
            logger.info("Creating application file {}".format(p))
            os.symlink(outfile.replace(".gz", ""), p)
    return

def install_annotation(build):
    """Install annotation files"""
    url = config[build]['annotation']['url']
    chr = config[build]['chr']
    builddir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', config[build]['species'], build)
    dest = os.path.join(builddir, "annotation", config[build]['annotation']['dest'])
    # Add annotation gtf 
    p = os.path.join(builddir, "annotation")
    if not os.path.exists(p):
        logger.info("Creating annotation directory {}".format(p))
        os.mkdir(p)
    p = os.path.join(builddir, "annotation", "gencode.v7.annotation.{}.gtf".format(chr))
    if not os.path.exists(p):
        logger.info("Creating annotation file {}".format(p))
        try:
            download = os.path.join(os.path.dirname(p), os.path.basename(url))
            logger.info("Downloading {} from {} with curl".format(download, url))
            cl = ["curl", url, "-o", download + ".tmp"]
            p = Popen(cl, stdout=PIPE)
            stdout, stderr = p.communicate()
            if p.returncode == 0:
               logger.info("{} generation completed; renaming to {}".format(download + ".tmp", download))
               os.rename(download + ".tmp", download)
            else:
                raise Exception("{} failed: \n{}".format(cmd, " ".join([stderr])))
            cmd = 'zcat {} | head -{} > {}'.format(download, config[build]['annotation']['lines'], dest)
            pipe = Popen(['/bin/bash', '-c', cmd])
            stdout, stderr = pipe.communicate()
            if pipe.returncode == 0:
                logger.info("{} generation completed; removing {}".format(dest, download))
                os.unlink(download)
            else:
                raise Exception("{} failed: \n{}".format(cmd, " ".join([stderr])))
        except:
            logger.info("Failed to download {}".format(download))


def setUp():
    """Set up test fixtures"""
    logger.info("Setting up test fixtures")
    install_ucsc_genome(build="hg19", start=config["hg19"]["start"], end=config["hg19"]["end"])
    install_annotation(build="hg19")
    
