# Copyright (C) 2015 by Per Unneberg
import re
import os
from snakemakelib.config import get_sml_config
from snakemakelib.utils import isoformat

sml_config = get_sml_config()

# SAM format specification
# @RG Read group. Unordered multiple @RG lines are allowed.

#   ID* Read group identifer. Each @RG line must have a unique ID. The
#   value of ID is used in the RG tags of alignment records. Must be
#   unique among all read groups in header section. Read group IDs may
#   be modied when merging SAM les in order to handle collisions.

# CN Name of sequencing center producing the read.
# DS Description.
# DT Date the run was produced (ISO8601 date or date/time).

# FO Flow order. The array of nucleotide bases that correspond to the
# nucleotides used for each ow of each read. Multi-base ows are
# encoded in IUPAC format, and non-nucleotide ows by various other
# characters. Format: /\*|[ACMGRSVTWYHKDBN]+/

# KS The array of nucleotide bases that correspond to the key sequence of each read.

# LB Library.

# PG Programs used for processing the read group.

# PI Predicted median insert size.

# PL Platform/technology used to produce the reads. Valid values:
# CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and
# PACBIO.

# PM Platform model. Free-form text providing further details of the platform/technology used.

# PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.

# SM Sample. Use pool name where a pool is being sequenced.

_read_group_dict =  {'ID':'identifier', 'CN':'center', 'DS':'description', 'DT':'date', 'FO':'floworder', 'KS':'keysequence', 'LB':'library', 'PG':'program', 'PI':'insertsize', 'PL': }
def read_group_from_str(s):
    """Generate read group values from a string.

    Tries to guess sensible values from string.

    Args: 
      s: <string> to be parsed

    Returns:
      d: <dict> of read group key:value mappings
    """
    ngs_cfg = get_sml_config('bio.ngs.settings')
    if not len(ngs_cfg['run_id_re']) == 2:
        raise IndexError("bio.ngs.settings.run_id_re key must be a tuple of length 2! Either reimplement the read group function or set the platform")
    m = re.match(ngs_cfg['run_id_re'][1], s)
    d = {k:"" for k in ngs_cfg['read_group_keys']}
    d.update({k:v for (k,v) in (zip(ngs_cfg['run_id_re'][0], m.groups())) if k in ngs_cfg['read_group_keys']})
    d['date'] = isoformat(d['date'])
    if "center" in list(d):
        d["center"] = ngs_cfg.get('center', "")
    if "platform" in list(d):
        d["platform"] = ngs_cfg.get('platform', "")
    if "description" in list(d):
        d["description"] = s
    if "id" in list(d):
        d["id"] = s
    return d

def find_files(path, re_str):
    """Find files in path that comply with a regular expression.

    Args:
      path:   path to search
      re_str: regular expression string

    Returns:
      flist: list of file names
    """
    r = re.compile(re_str)
    flist = []
    for root, dirs, files in os.walk(path):
        flist += [os.path.join(root, x) for x in files if r.match(x)]
    return sorted(flist)
