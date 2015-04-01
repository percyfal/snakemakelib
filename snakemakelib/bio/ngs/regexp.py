# Copyright (C) 2015 by Per Unneberg
import re
import os
from itertools import groupby
from snakemakelib.utils import isoformat
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

# Match beginning of path if it starts with dot or slash(?)
REGEXP_DOT_MATCH = r"(?:[\.\w\/]+)?\/"
# Match generic spacer sequence of >=0 characters
REGEXP_SPACER_MATCH = r"(?:.*)"

class MissingRequiredKeyException(Exception):
    """Exception if required key is missing"""

class FormatException(Exception):
    """Exception for malformatted entry"""

class UnimplementedException(Exception):
    """Exception for unimplemeted method"""

class DisallowedKeyException(Exception):
    """Exception for disallowed key"""

class RegexpDict(dict):
    _required_keys = []
    _group_keys = []
    _extra_keys = []

    def __init__(self, regexp=None, concat="_", *args, **kwargs):
        super(RegexpDict, self).__init__()
        # Set key values if in kwargs
        self.update({k:v  for k,v in kwargs.items() if k in self._group_keys})
        self._concat = concat
        self._init_regexp(regexp)
        self._init_format()
        
    def _init_regexp(self, regexp):
        self._regexp = re.compile(regexp)
        # Group keys according to prefix, e.g. PU=[PU1, PU2, PU3], and
        # SM = [SM]
        keymap = sorted([(re.sub("[0-9]+$", "", k), k)  if re.search("[0-9]+$", k) else (k, k) for k in list(self.re.groupindex.keys())])
        self._keymap = {k:[y[1] for y in list(v)] for (k,v) in groupby(keymap, lambda x: x[0])}
        self._validate_keys()

    @property
    def re(self):
        return self._regexp
        
    @property
    def pattern(self):
        return self.re.pattern

    @property
    def basename_pattern(self):
        """Return the basename pattern, replacing ?P= expressions with the corresponding group expression"""
        invmap = {v:k for k,v in self.re.groupindex.items()}
        remap = {k:"" for k,v in self.re.groupindex.items()}
        i = 1
        for m in re.finditer(r"\(?P<[A-Za-z0-9]+>([^\)]+)\)", self.pattern):
            remap[invmap[i]] = m.group(1)
            i += 1
        fmt = re.sub(r"\(?P=([A-Za-z0-9]+)\)", "P<\\1>{\\1})", self.pattern)
        return os.path.basename(fmt.format(**remap))

    @property
    def fmt(self):
        return self._fmt

    @property
    def allowed_keys(self):
        return list(set(self._required_keys + self._group_keys + self._extra_keys))
    
    def _init_format(self):
        """Initialize formatting string. Find groups defined by (?P<GROUP>) as
        well as constant expressions and concatenate"""
        m = re.findall("(\(\?P[<=](\w+)>?|({sep})|(?:[\[\]A-Za-z0-9\-\+\_]+\))|([A-Za-z0-9]+))".format(sep=os.sep), self.pattern)
        fmtlist = []
        for x in m:
            if x[1]:
                fmtlist.append("{" + x[1] + "}")
            elif x[2]:
                fmtlist.append(x[2])
            elif x[3]:
                fmtlist.append(x[3])
        self._fmt = re.sub("_{sep}_".format(sep=os.sep), os.sep, ("_".join(fmtlist)))



    def _validate_keys(self):
        """Validate keys. Importantly, make sure keys with indices, i.e. keys
        KEY1, KEY2, etc, stripped of indices are in the allowed keys.
        The sole purpose of the indices is to make the keys unique
        """
        seen_required = False
        for key, group in self._keymap.items():
            if not key in self.allowed_keys:
                raise DisallowedKeyException("key {key} not in allowed key set {allowed}".format(key=key, allowed=self.allowed_keys))
            if key in self._required_keys:
                seen_required = True
        if self._required_keys and not seen_required:
            raise MissingRequiredKeyException("one or several of the required keys '{keys}' is/are missing".format(keys=self._required_keys))

    def _parse_str(self, s, suffix):
        """Parse string and set read group dict"""
        pattern = self.pattern
        if s.startswith(os.curdir) or s.startswith(os.sep):
            pattern = REGEXP_DOT_MATCH + pattern
        if suffix:
            pattern = pattern + REGEXP_SPACER_MATCH + suffix
        m = re.match(pattern, s)
        if m is None:
            smllogger.debug("Unable to parse string {s} with regexp {re}".format(s=os.path.basename(s), re=self.re.pattern))
            return
        # Regular keys
        self.update({k:v for (k,v) in m.groupdict().items() if k not in self._extra_keys})
        self._process_indexed_keys(m)
        self._post_process_keys(m)

    def _process_indexed_keys(self, m):
        """Process indexed keys to unindexed version, e.g. collect PU1, PU2 to PU."""
        # Add indexed keys
        for key, group in self._keymap.items():
            if key in self._extra_keys:
                continue
            self.update({key: self._concat.join(m.group(mkey) for mkey in group)})

    def _post_process_keys(self, m):
        pass

    def parse(self, s, suffix=""):
        """Parse string and return string representation.

        Args:
          s (string): string to parse
          suffix (string): suffix to add to regular expression in search

        Returns:
          Regexp object
        """
        self.clear()
        self._parse_str(s, suffix)
        return self

class SampleRegexp(RegexpDict):
    _required_keys = ['SM']
    _group_keys = ['PU']
    _extra_keys = ['PATH']

    def __init__(self, regexp=None, *args, **kwargs):
        super(SampleRegexp, self).__init__(regexp, *args, **kwargs)

    def _post_process_keys(self, m):
        self['PATH'] = os.path.dirname(m.string)

class RunRegexp(RegexpDict):
    _group_keys = ['SM', 'PU', 'DT']
    _extra_keys = ['PATH']

    def __init__(self, regexp=None, *args, **kwargs):
        super(RunRegexp, self).__init__(regexp, *args, **kwargs)

    def _post_process_keys(self, m):
        self['PATH'] = os.path.dirname(m.string)

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

class ReadGroup(RunRegexp):
    """Adds formatting function for generating read group option string"""
    _group_keys = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
    _extra_keys = ['PATH']
    _group_dict =  {'ID' : 'identifier', 'CN' : 'center', 'DS' : 'description', 'DT' : 'date', 'FO' : 'floworder', 'KS' : 'keysequence', 'LB' : 'library', 'PG' : 'program', 'PI' : 'insertsize', 'PL': 'platform', 'PU' : 'platform-unit', 'SM' : 'sample'}

    def __init__(self, regexp=None, opt_prefix="--", *args, **kwargs):
        super(ReadGroup, self).__init__(regexp, *args, **kwargs)
        self._opt_prefix = opt_prefix

    def _post_process_keys(self, m):
        self['PATH'] = os.path.dirname(m.string)
        if not 'ID' in self.keys() or not self.get('ID', ""):
            inv_map = {v:k for (k,v) in list(self.re.groupindex.items())}
            self['ID'] = os.path.basename(self.fmt.format(**self))

    def _fmt_string(self, k):
        """Take care of date string"""
        if k == 'DT':
            return isoformat(self[k])
        return self[k]
        
    def __str__(self):
        """Return a generic program string"""
        return " ".join(["{dash}{key} {value}".format(dash=self._opt_prefix, key=self._group_dict[k], value=self._fmt_string(k)) for k in sorted(list(self.keys())) if not self[k] is None and k in self._group_keys])

    
