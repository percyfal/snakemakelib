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
        import time
        logger.warn ("\n\nno metadata file '{metadata}' found;\n\nplease initiate analysis by running 'snakemake {metadata}'\n\n".format(metadata=metadata_file))
        time.sleep(3)
    return metadata_list
