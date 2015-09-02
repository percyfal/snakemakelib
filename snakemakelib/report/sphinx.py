# Copyright (c) 2014 Per Unneberg
import os
import re
import glob
import texttable as tt
from collections import namedtuple
from mako.template import Template
from snakemakelib.report.picard import PicardMetricsCollection
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

TEMPLATEPATH = os.path.join(os.path.dirname(__file__), os.pardir, "data", "templates", "doc")

templates = {
    'make' : (Template(filename=os.path.join(TEMPLATEPATH, "Makefile.mako")),''),
    'sample' : (Template(filename=os.path.join(TEMPLATEPATH, "source", "samples", "sample.mako")), '.rst'),
    'index' : (Template(filename=os.path.join(TEMPLATEPATH, "source", "index.mako")), '.rst'),
#    'sampleindex' : (Template(filename=os.path.join(TEMPLATEPATH, "source", "samples", "index.mako")), '.rst'),
    'conf' : (Template(filename=os.path.join(TEMPLATEPATH, "source", "conf.mako")), '.py')
    }

docdirs = ["build", "source", os.path.join("source", "_static"), os.path.join("source", "_templates"),
           os.path.join("source", "samples")]

def _setup_directories(root):
    """create directory structure"""
    if not os.path.exists(root):
        smllogger.info("Making documentation output directory '{}'".format(root))
        for d in docdirs:
            if not os.path.exists(os.path.join(root, d)):
                smllogger.info("Making directory " + str(os.path.join(root, d)))
                os.makedirs(os.path.join(root, d))

def make_rst_table(data):
    """Make rst table with :py:mod:`Texttable`"""
    if data is None:
        return ""
    else:
        tab_tt = tt.Texttable()
        tab_tt.set_precision(2)
        tab_tt.add_rows(data)
        return tab_tt.draw()


def sphinx_sample_metrics_report(input, config):
    Sample = namedtuple('Sample', ['sample_id', 'project_id', 'pmc'])
    metrics = {}
    for s in input:
        mlist = [(os.path.basename(s), m) for m in glob.glob(os.path.join(s, "*metrics"))]
        pmc = PicardMetricsCollection(mlist)
        st = Sample(sample_id=os.path.basename(s), project_id=config['project_id'], pmc=pmc)
        metrics[st.sample_id] = st
    # Setup output directory
    kw = {
        'docroot':config['docroot'],
        'sample_id':None,
        'metrics':{k: (st.sample_id, st.project_id, st.pmc.metrics(as_csv=True)) for k,st in metrics.items()},
        'project_id':config['project_id'],
        'project_name':config['project_name'],
        'application':config['application'],
        'date':config['date'],
        'samples': "\n".join(["   {}".format(s.sample_id) for s in sorted(metrics.values())])
    }
    _setup_directories(kw['docroot'])
    for k,v in templates.items():
        if k == "sample":
            continue
        else:
            outfile = os.path.splitext(os.path.join(kw['docroot'], os.path.relpath(v[0].filename, TEMPLATEPATH)))[0] + v[1]
            with open(outfile, "w") as fh:
                fh.write(v[0].render(**kw))
    for s in sorted(metrics.values()):
        kw['sample_id'] = s.sample_id
        kw['pmcmetrics'] = s.pmc.metrics(as_csv=True)
        kw['pmcidlist'] = s.pmc.idlist()
        outfile = os.path.join(kw['docroot'], os.path.dirname(os.path.relpath(templates['sample'][0].filename, TEMPLATEPATH)), "{}.rst".format(s.sample_id))
        with open(outfile, "w") as fh:
            fh.write(templates['sample'][0].render(**kw))
