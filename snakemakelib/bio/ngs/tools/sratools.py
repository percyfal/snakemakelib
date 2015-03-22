# Copyright (C) 2015 by Per Unneberg

import os
import csv
from snakemakelib.log import LoggerManager

logger = LoggerManager().getLogger(__name__)

def get_metadata_list(metadata_file):
    """
    Read an SRA project file and return entries as a list. Will issue a
    warning if file does not exists.

    Args: 
      metadata - file name

    Returns:
      list of file contents where every item is a dictionary

    """
    metadata_list = []
    if os.path.exists(metadata_file):
        with open(metadata_file, "r") as fh:
            reader = csv.DictReader(fh.readlines())
        metadata_list = [row for row in reader]
    else:
        import sys
        if metadata_file in sys.argv:
            return metadata_list
        import time
        logger.warn("""

        no metadata file '{metadata}' found

        please initiate analysis by running 'snakemake {metadata}'

        """.format(metadata=metadata_file))
        time.sleep(3)
    return metadata_list
