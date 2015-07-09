# Copyright (C) 2015 by Per Unneberg
import re
import os
import csv
from snakemakelib.bio.ngs.regexp import RegexpDict
from snakemakelib.bio.ngs.utils import find_files
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)


def generic_target_generator(tgt_re, src_re=None, samples=[], runs=[],
                             sample_column_map={}, sampleinfo="",
                             target_suffix="", filter_suffix="", **kwargs):
    """Generic target generator.

    Args:
      tgt_re (RegexpDict): RegexpDict object corresponding to the target
                           regular expression
      src_re (RegexpDict): RegexpDict object corresponding to the source
                           regular expression
      samples (list): sample names
      runs (list): run names
      sample_column_map (dict): mapping from sampleinfo column names to
                         regexp group names, e.g.
                         {'SampleID':'SM', 'Lane':'PU1'}
      sampleinfo (str): sample information file
      target_suffix (str): suffix of generated targets
      filter_suffix (str): suffix to use for filtering when generating target
                     names based on input files

    Returns:
      list of target names
    """
    assert isinstance(tgt_re, RegexpDict),\
        "tgt_re argument must be of type {}".format(RegexpDict)
    if src_re is None:
        src_re = tgt_re
    assert isinstance(src_re, RegexpDict),\
        "src_re argument must be of type {}".format(RegexpDict)
    # 1. Generate targets from command line options
    if samples and runs:
        smllogger.debug("trying to gather target information based on " +
                        "configuration keys 'samples' and 'runs'")
        if len(samples) == len(runs):
            cfg_list = list(zip(samples, runs))
            mlist = []
            for (s, r) in cfg_list:
                # Use basename searches for samples and runs
                m = re.search(src_re.basename_pattern, r).groupdict()\
                    if not re.search(src_re.basename_pattern, r) is None else {}
                if m:
                    m.update({'SM': s})
                    mlist.append(m)
            tgts = [tgt_re.fmt.format(**ml) + target_suffix for ml in mlist]
            return sorted(tgts)
        else:
            smllogger.warn("if samples and runs are provided, they must be of equal lengths")

    # 2. Generate targets from information in samplesheet
    if sampleinfo != "":
        smllogger.debug("trying to gather target information from configuration key 'sampleinfo'")
        if isinstance(sampleinfo, str) and not os.path.exists(sampleinfo):
            smllogger.debug("no such sample information file '{sampleinfo}'; trying to deduct targets from existing files".format(sampleinfo=sampleinfo))
        else:
            smllogger.debug("Reading sample information from '{sampleinfo}'".format(sampleinfo=sampleinfo))
            if isinstance(sampleinfo, str):
                with open(sampleinfo, 'r') as fh:
                    reader = csv.DictReader(fh.readlines())
            else:
                reader = sampleinfo
                assert type(reader) is csv.DictReader,\
                    "sampleinfo is not a 'csv.DictReader'; if not a file name, must be a 'csv.DictReader'"
            reader.fieldnames = [fn if fn not in sample_column_map.keys()
                                 else sample_column_map[fn]
                                 for fn in reader.fieldnames]
            if samples:
                tgts = [tgt_re.fmt.format(**row) + target_suffix
                        for row in reader if row['SM'] in samples]
            else:
                tgts = [tgt_re.fmt.format(**row) + target_suffix
                        for row in reader]
            return sorted(tgts)

    # 3. Generate targets from input files
    smllogger.debug("Getting sample information from input files")
    inputs = find_files(regexp=src_re.basename_pattern + filter_suffix,
                        limit={'SM': samples} if samples else {})
    if inputs:
        tgts = [tgt_re.fmt.format(**src_re.parse(f))
                + target_suffix for f in inputs]
        return sorted(tgts)
    smllogger.warn("No targets could be generated!")
    return []
