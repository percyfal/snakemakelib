# Copyright (c) 2014 Per Unneberg
import os
import glob
import itertools
from collections import OrderedDict
from snakemakelib.report.picard import PicardMetricsCollection

def group_samples(samples, grouping="sample"):
    """Group samples by sample or sample run.

    Args:
      samples: list of namedtuples with fields (sample_id, project_id)
      grouping: what to group by

    Returns:
      groups: dictionary of grouped items
    """
    groups = OrderedDict()
    if grouping == "sample":
        for k,g in itertools.groupby(samples, key=lambda x:x.sample_id()):
            groups[k] = list(g)
    return groups
