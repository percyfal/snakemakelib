# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import unittest
import logging
from snakemakelib.config import sml_config
from collections import OrderedDict, namedtuple

logging.basicConfig(level=logging.DEBUG)

config = dict()

class BaseConfig(dict):
    _sections = []
    def __init__(self, *args, **kwargs):
        self._sections = [kk for k in args for kk in list(k)] + list(kwargs)
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        if not key in self._sections:
            raise KeyError("key '" + key + "' not found in configuration dictionary")
        dict.__setitem__(self, key, val)

d = {'test':1, 't2': 2}
b = BaseConfig(d, a=2)
print (b)
#b['23'] = "test"
#print (b)

# print (b._sections)
# print (dir(b._sections))
# print (type(b._sections))
# print ('test' in b._sections)
# for k in b._sections:
#     print (k)

def fun1():
    return "fun1"

def fun2():
    return "fun2"

def fun3():
    return " ".join(["-R", str(cfg['bio.ngs.settings']['db']['ref'])])

c = BaseConfig({'bio.ngs.settings' : 
                BaseConfig({'fastq_suffix':None,
                            'flowcells': None,
                            'flowcellruns' : None,
                            'annotation': BaseConfig({'annot_label':None, 'transcript_annot_gtf':None}),
                            'db' : BaseConfig({'dbsnp':None, 'ref':None})})
                })
c['bio.ngs.settings']['fastq_suffix'] = ".fastq.gz"
c['bio.ngs.settings']['db']['ref'] = "hg19"
d = BaseConfig({'bio.ngs.align.bwa' : 
                BaseConfig({'cmd':'bwa',
                            'ref': None,
                            'threads' : None,
                            'options': "-M",
                            'mem' : fun1(),
                        })
                })
global cfg
cfg = BaseConfig({'bio.ngs.settings':c['bio.ngs.settings'],
                         'bio.ngs.align.bwa':d['bio.ngs.align.bwa']})


class TestBasicConfig(unittest.TestCase):
    def test_print_config(self):
        """Test printing config"""
        print ("Configuration", config)

    def test_sml_config(self):
        """Test printing global sml config"""
        print ("Global configuration " , sml_config)

    def test_config(self):
        print (c)


    def test_base_config(self):

        print (cfg)
        print (cfg['bio.ngs.settings'])
        print (cfg['bio.ngs.align.bwa'])
        print (cfg['bio.ngs.align.bwa']['mem'])
        cfg['bio.ngs.align.bwa']['mem'] = fun3()
        print (cfg['bio.ngs.align.bwa']['mem'])
