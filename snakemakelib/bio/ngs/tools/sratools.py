# Copyright (C) 2015 by Per Unneberg
import os
import csv
from snakemakelib.config import update_snakemake_config
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def register_metadata(metadata_file, config):
    """Read an SRA project file and register metadata in sml_config. Will
    issue a warning if file does not exists.

    Args: 
      metadata - file name
      config - configuration to update
    """
    metadata_list = []
    import sys
    if metadata_file in sys.argv:
        return metadata_list
    try:
        with open(metadata_file, "r") as fh:
            reader = csv.DictReader(fh.readlines())
        metadata_list = [row for row in reader]
        run2sample = {row["Run"]:row["SampleName"] for row in metadata_list}
        config = update_snakemake_config(config, {
            'bio.ngs.settings' : {'sampleinfo' : metadata_file},
            'bio.ngs.tools.sratools': {'_datadir': os.path.dirname(metadata_file),
                                       '_run2sample' : run2sample,
                                       '_metadata' : metadata_list}})
    except Exception:
        raise Exception("""

        no metadata file '{metadata}' found

        please initiate analysis by running 'snakemake {metadata}'

        """.format(metadata=metadata_file))

    return config
