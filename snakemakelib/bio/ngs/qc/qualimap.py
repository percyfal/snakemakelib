# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import jinja2
import numpy as np
from snakemakelib.results import Results
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

class Qualimap(Results):
    
    def __init__(self, *args, **kw):
        super(Qualimap, self).__init__()

    def _collect_results(self):
        smllogger.info("Collecting results")
        for (f, s) in zip(self._inputfiles, self._samples):
            df_tmp = 

def collect_qualimap_results(inputfiles, samples):
    """Collect qualimap results.
    
    Args:
      inputfiles (list): list of input file names
      samples (list): list of sample names

    Returns:
      dict: dictionary of pandas data frames

    """
    pass
