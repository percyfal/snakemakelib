# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import gzip
import io
import unittest
import logging
import shutil
import subprocess
from nose.tools import raises
from nose.plugins.attrib import attr

logger = logging.getLogger(__name__)

BASEPATH = os.path.dirname(__file__)
ROOTPATH = os.path.join(BASEPATH, os.pardir, os.pardir, os.pardir)

snakefile = """# -*- snakemake -*-
import os
import sys
import re
from snakemakelib.config import update_snakemake_config, load_sml_config
from snakemakelib.utils import rreplace
from snakemakelib.bio.ngs.targets import generic_target_generator

workdir: '{workdir}'

config = load_sml_config(config)

local_config = {{
    'bio.ngs.settings' : {{
        'sample_organization' : 'Illumina@SciLife',
        'samples' : config.get("samples", ['P001_101', 'P001_102']),
        'runs' : config.get("runs", []),
        'db' : {{
            'ref' : '{ref}',
            'dbsnp' : '{dbsnp}',
        }},
        'java' : {{
            'java_mem' : '1g',
        }},
        'sequence_capture' : {{
            'bait_regions' : '{baits}',
            'target_regions' : '{targets}',
        }},
    }},
    'bio.ngs.align.bwa' : {{
        'index' : '{bwaref}',
    }},
    'bio.ngs.qc.picard' : {{
        'home' : '{picard}',
    }},
    'bio.ngs.tools.gatk' : {{
        'home' : '{gatk}',
    }},
    'bio.ngs.variation.variation' : {{
        'snpeff' : {{
             'path' : '{snpeff}',
             'genome_version' : 'hg19',
        }},
    }},
}}

config = update_snakemake_config(config, local_config)

include: '{variation}'

cfg = config['bio.ngs.settings']
path = cfg.get('path') if not cfg.get('path') is None else os.curdir

"""
         
@unittest.skipIf((os.getenv("PICARD_HOME") is None or os.getenv("PICARD_HOME") == ""), "No Environment PICARD_HOME set; skipping")
@unittest.skipIf((os.getenv("GATK_HOME") is None or os.getenv("GATK_HOME") == ""), "No Environment GATK_HOME set; skipping")
@unittest.skipIf((os.getenv("SNPEFF_HOME") is None or os.getenv("SNPEFF_HOME") == ""), "No Environment SNPEFF_HOME set; skipping")
@unittest.skipIf(shutil.which('bwa') is None, "No executable bwa found; skipping")
@unittest.skipIf(shutil.which('samtools') is None, "No executable samtools found; skipping")

def setUp():
    logger.info("Setting up text fixtures for {}".format(__name__))
    
    # Create Snakefile in test directory
    with open("Snakefile", "w") as fh:
        fh.write(snakefile.format(workdir=os.path.join(ROOTPATH, 'data', 'projects', 'J.Doe_00_01'),
                                  variation=os.path.join(ROOTPATH, 'snakemakelib', 'workflows', 'bio', 'variation.workflow'),
                                  ref=os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'seq', 'chr11.fa'),
                                  dbsnp=os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'variation', 'dbsnp132_chr11.vcf'),
                                  bwaref=os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'bwa', 'chr11.fa'),
                                  picard=os.getenv("PICARD_HOME"),
                                  gatk=os.getenv("GATK_HOME"),
                                  snpeff=os.getenv("SNPEFF_HOME"),
                                  baits=os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'seqcap', 'chr11_baits.interval_list'),
                                  targets=os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'seqcap', 'chr11_targets.interval_list')
                              ))
    # Format bwa ref if not already done
    genomedir = os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'seq')
    bwadir = os.path.join(ROOTPATH,  'data', 'genomes', 'Hsapiens', 'hg19', 'bwa')
    if not os.path.exists(os.path.join(bwadir, 'chr11.fa.bwt')):
        subprocess.check_call(['bwa', 'index', os.path.join(bwadir, 'chr11.fa')])
    if not os.path.exists(os.path.join(genomedir, 'chr11.dict')):
        subprocess.check_call(['java', '-jar', os.path.join(os.getenv("PICARD_HOME"), 'picard.jar'), 'CreateSequenceDictionary', "R={}".format(os.path.join(genomedir, 'chr11.fa')), "O={}".format(os.path.join(genomedir, 'chr11.dict'))])
    if not os.path.exists(os.path.join(genomedir, 'chr11.fa.fai')):
        subprocess.check_call(['samtools', 'faidx', os.path.join(genomedir, 'chr11.fa')])

@attr('slow')
class TestVariationWorkflow(unittest.TestCase):
    def test_variation_workflow(self):
        """Test variation.workflow.

        Don't run snpEff as this step may fail on low-memory machines,
        such as laptops.
        """
        outputs = ['P001_101/P001_101.sort.merge.rg.dup.realign.recal.bp_variants.phased.vcf',
                   'P001_102/P001_102.sort.merge.rg.dup.realign.recal.bp_variants.phased.vcf']
        subprocess.check_call(['snakemake', '-F'] + outputs)
        subprocess.check_call(['snakemake', '-F', 'variation_metrics'])

class TestBwaAlign(unittest.TestCase):
    def test_bwa_align(self):
        """Test bwa alignment"""
        bwaout = 'P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam'
        subprocess.check_call(['snakemake', '-F', bwaout])
        subprocess.check_call(['rm', '-f', os.path.join(ROOTPATH,  'data', 'projects', 'J.Doe_00_01', 'P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam')])
