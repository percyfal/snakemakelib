# Copyright (C) 2015 by Per Unneberg
import re
import os
import csv
from snakemakelib.bio.ngs.utils import find_files
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def generic_target_generator(fmt, rg, cfg, path=os.curdir, prepend_path=True):
    """Generic target generator.

    Args:

      fmt: python miniformat string detailing what the target should
           look like. The format names are based on the ReadGroup
           identifiers. Example: "{SM}/{PU}/{PU}_{SM}_1.fastq.gz will
           generate a target residing in path SM, with subdirectory PU
           (platform unit), and named platform unit underscore sample
           underscore .fastq.gz.

      rg: ReadGroup object specifying how format names are derived
          from string

      cfg: Configuration dictionary for bio.ngs.settings

      path: path to search in; usually the snakemake workdir

      prepend_path: prepend path to the targets

    Returns:
      targets: list of target names

    """
    if prepend_path:
        ppath = path if path != os.curdir else ""
    else:
        ppath = ""
    # 1. from command line options
    if cfg['samples'] and cfg['runs']:
        if not len(cfg['samples']) == len(cfg['runs']):
            raise Exception("if samples and runs are provided, they must be of equal lengths")
        cfg_list = list(zip(cfg['samples'], cfg['runs']))
        mlist = []
        for (s, r) in cfg_list:
            m = re.match(rg.pattern, r).groupdict() if not re.match(rg.pattern, r) is None else {}
            if m:
                m.update({'SM':s})
                mlist.append(m)
        tgts = [fmt.format(**m) for m in mlist]
        return [os.path.join(ppath, t) for t in tgts]

    # 2. Read samplesheet here
    if cfg['sampleinfo'] != "":
        if isinstance(cfg['sampleinfo'], str) and not os.path.exists(cfg['sampleinfo']):
            smllogger.info("no such sample information file '{sampleinfo}'; trying to deduct targets from existing files".format(sampleinfo=cfg['sampleinfo']))
        else:
            smllogger.info("Reading sample information from '{sampleinfo}'".format(sampleinfo=cfg['sampleinfo']))
            if isinstance(cfg['sampleinfo'], str):
                with open(cfg['sampleinfo'], 'r') as fh:
                    reader = csv.DictReader(fh.readlines())
            else:
                reader = cfg['sampleinfo']
                assert type(reader) is csv.DictReader, "cfg['sampleinfo'] is not a 'csv.DictReader'"
            reader.fieldnames = [fn if fn != cfg['sample_column_name'] else 'SM' for fn in reader.fieldnames]
            reader.fieldnames = [fn if fn != cfg['run_column_name'] else 'PU' for fn in reader.fieldnames]
            if cfg['samples']:
                tgts = [fmt.format(**row) for row in reader if row['SM'] in cfg['samples']]
            else:
                tgts = [fmt.format(**row) for row in reader]
            return [os.path.join(ppath, t) for t in tgts]

    # 3. generate from input files
    limit = {}
    if cfg['samples']:
        limit['SM'] = cfg['samples']
    inputs = find_files(path=path, re_str=rg.pattern, limit=limit)
    if inputs:
        rgfmt = [dict(rg.parse(f)) for f in inputs]
        tgts = [fmt.format(**f) for f in rgfmt]
        return [os.path.join(ppath, t) for t in tgts]
    return []

    
