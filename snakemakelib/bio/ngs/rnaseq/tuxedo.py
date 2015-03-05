# Copyright (C) 2015 by Per Unneberg
from snakemakelib.config import get_sml_config
from snakemakelib.bio.ngs.utils import read_group_from_str

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
def opt_read_group(prefix):
    """Generate read group string for tuxedo alignments.

    Tries to guess sensible values from the file prefix.

    Args:
      prefix: <string> to be parsed

    Returns:
      opt: tophat-formatted <string> to be used in options
    """
    d = read_group_from_str(prefix)
    s = " ".join(["--rg-{k} {v}".format(k=k, v=v) for (k,v) in sorted(d.items()) if not v == ""])
    return s
    
