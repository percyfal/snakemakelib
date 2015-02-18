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

logger = logging.getLogger(__name__)

snakefile = """# -*- snakemake -*-
import os
from snakemakelib.config import init_sml_config, get_sml_config

workdir: '{workdir}'

local_config = {{
    'bio.ngs.settings' : {{
        'fastq_suffix' : ".fastq.gz",
        'samples' : config.get("samples", ['P001_101', 'P001_102']),
        'flowcells' : config.get("flowcells", ['120924_AC003CCCXX', '121015_BB002BBBXX']),
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
        'ref' : '{bwaref}',
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

init_sml_config(local_config)

include: '{variation}'

"""
         
@unittest.skipIf((os.getenv("PICARD_HOME") is None or os.getenv("PICARD_HOME") == ""), "No Environment PICARD_HOME set; skipping")
@unittest.skipIf((os.getenv("GATK_HOME") is None or os.getenv("GATK_HOME") == ""), "No Environment GATK_HOME set; skipping")
@unittest.skipIf((os.getenv("SNPEFF_HOME") is None or os.getenv("SNPEFF_HOME") == ""), "No Environment SNPEFF_HOME set; skipping")
@unittest.skipIf(shutil.which('bwa') is None, "No executable bwa found; skipping")
#@unittest.skipIf(shutil.which('bowtie') is None, "No executable bowtie found; skipping")
@unittest.skipIf(shutil.which('samtools') is None, "No executable samtools found; skipping")

def setUp():
    logger.info("Setting up text fixtures for {}".format(__name__))
    
    # Create Snakefile in test directory
    with open("Snakefile", "w") as fh:
        fh.write(snakefile.format(workdir=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'projects', 'J.Doe_00_01'),
                                  variation=os.path.join(os.path.abspath(os.curdir), os.pardir, 'workflows', 'bio', 'variation.workflow'),
                                  ref=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'seq', 'chr11.fa'),
                                  dbsnp=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'variation', 'dbsnp132_chr11.vcf'),
                                  bwaref=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'bwa', 'chr11.fa'),
                                  picard=os.getenv("PICARD_HOME"),
                                  gatk=os.getenv("GATK_HOME"),
                                  snpeff=os.getenv("SNPEFF_HOME"),
                                  baits=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'seqcap', 'chr11_baits.interval_list'),
                                  targets=os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'seqcap', 'chr11_targets.interval_list')
                              ))
    # Format bwa ref if not already done
    genomedir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'seq')
    bwadir = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'genomes', 'Hsapiens', 'hg19', 'bwa')
    if not os.path.exists(os.path.join(bwadir, 'chr11.fa.bwt')):
        subprocess.check_call(['bwa', 'index', os.path.join(bwadir, 'chr11.fa')])
    if not os.path.exists(os.path.join(genomedir, 'chr11.dict')):
        subprocess.check_call(['java', '-jar', os.path.join(os.getenv("PICARD_HOME"), 'picard.jar'), 'CreateSequenceDictionary', "R={}".format(os.path.join(genomedir, 'chr11.fa')), "O={}".format(os.path.join(genomedir, 'chr11.dict'))])
    if not os.path.exists(os.path.join(genomedir, 'chr11.fa.fai')):
        subprocess.check_call(['samtools', 'faidx', os.path.join(genomedir, 'chr11.fa')])

class TestVariationWorkflow(unittest.TestCase):
    def test_variation_workflow(self):
        """Test variation.workflow.

        Don't run snpEff as this step may fail on low-memory machines,
        such as laptops.
        """
        outputs = ['P001_101/P001_101.sort.merge.rg.dup.realign.recal.bp_variants.phased.vcf',
                   'P001_102/P001_102.sort.merge.rg.dup.realign.recal.bp_variants.phased.vcf']
        subprocess.check_call(['snakemake'] + outputs)
        subprocess.check_call(['snakemake', 'metrics'])


class TestBwaAlign(unittest.TestCase):
    def test_bwa_align(self):
        """Test bwa alignment"""
        bwaout = 'P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam'
        subprocess.check_call(['snakemake', '-F', bwaout])
        subprocess.check_call(['rm', '-f', os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'projects', 'J.Doe_00_01', 'P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam')])
