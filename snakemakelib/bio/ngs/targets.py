# Copyright (C) 2015 by Per Unneberg
import re
import os
import csv
from snakemakelib.bio.ngs.utils import find_files

def generic_target_generator(fmt, rg, cfg, path=os.curdir):
    """Generic target generator.

    Args:

      fmt: python miniformat string detailing what the target should
           look like. The format names are based on the ReadGroup
           identifiers. Example:
           "{PATH}/{SM}/{PU}/{PU}_{SM}_1.fastq.gz will generate a
           target residing in path PATH, with subdirectory SM (sample
           name), with subdirectory PU (platform unit), and named
           platform unit underscore sample underscore .fastq.gz.

      rg: ReadGroup object specifying how format names are derived
          from string

      cfg: Configuration dictionary for bio.ngs.settings

    Returns:
      targets: list of target names
    """
    # 1. from command line options
    if cfg['samples'] and cfg['flowcells'] and cfg['lanes']:
        if not len(cfg['samples']) == len(cfg['flowcells'] or len(cfg['samples'] == len(cfg['flowcells']))):
            raise Exception("if samples, flowcells, lanes all provided, must be of equal lengths")
        cfg_list = list(zip(cfg['samples'], cfg['flowcells'], cfg['lanes']))
        tgts = [fmt.format(SM=s, PATH=path, **dict(cfg)['platform_unit_fn']((s,fc,l))) for (s, fc, l) in cfg_list]
        return tgts

    # 2. Read samplesheet here
    if cfg['sampleinfo'] != "":
        if isinstance(cfg['sampleinfo'], str):
            with open(cfg['sampleinfo'], 'r') as fh:
                reader = csv.DictReader(fh.readlines())
        else:
            # Assume cvs.DictReader a
            reader = cfg['sampleinfo']
        assert type(reader) is csv.DictReader, "cfg['sampleinfo'] is not a 'csv.DictReader'"
        if cfg['samples']:
            tgts = [fmt.format(PATH=path, **row) for row in reader if row['SM'] in cfg['samples']]
        else:
            tgts = [fmt.format(PATH=path, **row) for row in reader]
        return tgts

    # 3. generate from input files
    inputs = find_files(path=path, re_str=rg.pattern)
    if inputs:
        tgts = [fmt.format(PATH=path, **dict(rg.parse(f))) for f in inputs]
        return tgts
    return []

    
