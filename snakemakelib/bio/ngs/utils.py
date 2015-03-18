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
#   be modified when merging SAMfiles in order to handle collisions.

# CN Name of sequencing center producing the read.
# DS Description.
# DT Date the run was produced (ISO8601 date or date/time).

# FO Flow order. The array of nucleotide bases that correspond to the
# nucleotides used for each ow of each read. Multi-base rows are
# encoded in IUPAC format, and non-nucleotide rows by various other
# characters. Format: /\*|[ACMGRSVTWYHKDBN]+/

# KS The array of nucleotide bases that correspond to the key sequence
# of each read.

# LB Library.

# PG Programs used for processing the read group.

# PI Predicted median insert size.

# PL Platform/technology used to produce the reads. Valid values:
# CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and
# PACBIO.

# PM Platform model. Free-form text providing further details of the platform/technology used.

# PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.

# SM Sample. Use pool name where a pool is being sequenced.

class MissingIDException(Exception):
    """Exception if read group ID is missing"""

class FormatException(Exception):
    """Exception for malformatted entry"""

class UnimplementedException(Exception):
    """Exception for unimplemeted method"""

class DisallowedKeyException(Exception):
    """Exception for disallowed key"""

class ReadGroup(dict):
    """Create a read group representation from string"""
    _read_group_keys = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
    _read_group_dict =  {'ID' : 'identifier', 'CN' : 'center', 'DS' : 'description', 'DT' : 'date', 'FO' : 'floworder', 'KS' : 'keysequence', 'LB' : 'library', 'PG' : 'program', 'PI' : 'insertsize', 'PL': 'platform', 'PU' : 'platform-unit', 'SM' : 'sample'}

    _extra_keys = ['PATH', 'NA']

    _allowed_keys = _read_group_keys + _extra_keys

    def __init__(self, run_id_re, opt_prefix="--", concat="_", *args, **kwargs):
        dict.__init__(self)
        self._run_id_re = run_id_re
        self._init_regex()
        self._opt_prefix = opt_prefix
        self._concat = concat
        self.update({x:None for x in self._read_group_keys})
        self.update(*args, **kwargs)

    def _init_regex(self):
        self._regex = re.compile(self._run_id_re)
        self._indexkeys = sorted([(re.sub("[0-9]+$", "", k), k) for k in list(self._regex.groupindex.keys()) if re.search("[0-9]+$", k)])
        # Validate group keys
        for k in self._regex.groupindex.keys():
            if not k in self._allowed_keys and not k in [v for (k,v) in self._indexkeys]:
                raise DisallowedKeyException("key {key} not in allowed key set {allowed}".format(key=k, allowed=self._allowed_keys))

        self._indexdict = {k:[] for (k,v) in self._indexkeys}
        [self._indexdict[k].append(v) for (k,v) in self._indexkeys]

    def _parse_str(self, s):
        """Parse string and set read group dict"""
        m = self._regex.match(os.path.basename(s))
        path = os.path.dirname(s)
        [self.update({k: m.group(k)}) for k in m.groupdict().keys() if not k in self._extra_keys]
        if self._indexdict:
            [self.update({k: self._concat.join(m.group(mkey) for mkey in self._indexdict[k])}) for k in self._indexdict.keys()]
        if not self['ID']:
            inv_map = {v:k for (k,v) in list(self._regex.groupindex.items())}
            self['ID'] = self._concat.join(m.group(i) for i in range(1, self._regex.groups + 1) if not inv_map[i] in self._extra_keys)
        self['PATH'] = path
            
    def parse(self, s):
        """Parse string and return string representation"""
        self._parse_str(s)
        return self
    
    @property
    def pattern(self):
        return self._regex.pattern

    def _validate_keys(self):
        # ID required!
        if self['ID'] is None:
            raise MissingIDException("Read group ID required")
        if not self['FO'] is None:
            if not re.search("\*|[ACMGRSVTWYHKDBN]+", self['FO']):
                raise FormatException("FO must be of format '\*|[ACMGRSVTWYHKDBN]+'")

    def _fmt(self, k):
        """Take care of date string"""
        if k == 'DT':
            return isoformat(self[k])
        return self[k]

    def __str__(self):
        """Return a generic program string"""
        self._validate_keys()
        return " ".join(["{dash}{key} {value}".format(dash=self._opt_prefix, key=self._read_group_dict[k], value=self._fmt(k)) for k in sorted(list(self.keys())) if not self[k] is None and k in self._read_group_keys])

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

def dict_to_R(d):
    """Translate a python dictionary to an R option string

    Args:
      d: dictionary

    Returns:
      A string representing an option string of the dictionary in R
    """
    outlist = []
    for (k,v) in d.items():
        if isinstance(v, list):
            if isinstance(v[0], int):
                outlist.append("{k} = c(".format(k=k) + ",".join("{vv}".format(vv=vv) for vv in v) + ")")
            else:
                outlist.append("{k} = c(".format(k=k) + ",".join("'{vv}'".format(vv=vv) for vv in v) + ")")     
        elif isinstance(v, dict):
            outlist.append("{k} = c(".format(k=k) + ",".join("{kk}={vv}".format(kk=kk, vv=vv) for (kk,vv) in v.items()) + ")")
        elif v is True:
            outlist.append("{k}=TRUE".format(k=k))
        elif v is False:
            outlist.append("{k}=FALSE".format(k=k))
        elif v is None:
            outlist.append("{k}=NULL".format(k=k))
        elif isinstance(v, int):
            outlist.append("{k}={v}".format(k=k, v=v))
        elif isinstance(v, float):
            outlist.append("{k}={v}".format(k=k, v=v))
        else:
            outlist.append("{k}='{v}'".format(k=k, v=v))
    return ",".join(outlist)
