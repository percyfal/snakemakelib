# Copyright (C) 2015 by Per Unneberg
from snakemakelib.config import get_sml_config
from snakemakelib.bio.ngs.utils import ReadGroup

sml_config = get_sml_config()

# From tophat2:
# --rg-id                        <string>    (read group ID)
# --rg-sample                    <string>    (sample ID)
# --rg-library                   <string>    (library ID)
# --rg-description               <string>    (descriptive string, no tabs allowed)
# --rg-platform-unit             <string>    (e.g Illumina lane ID)
# --rg-center                    <string>    (sequencing center name)
# --rg-date                      <string>    (ISO 8601 date of the sequencing run)
# --rg-platform                  <string>    (Sequencing platform descriptor)

class TuxedoReadGroup(ReadGroup):
    _read_group_dict =  {'ID' : 'id', 'CN' : 'center', 'DS' : 'description', 'DT' : 'date', 'FO' : 'floworder', 'KS' : 'keysequence', 'LB' : 'library', 'PG' : 'program', 'PI' : 'insertsize', 'PL': 'platform', 'PU' : 'platform-unit', 'SM' : 'sample'}

    def __init__(self, s, opt_prefix="--rg-", *args, **kwargs):
        ReadGroup.__init__(self, s, opt_prefix, *args, **kwargs)

