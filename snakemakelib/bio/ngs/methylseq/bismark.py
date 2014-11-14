# Copyright (C) 2014 by Per Unneberg
from snakemakelib.config import get_sml_config

sml_config = get_sml_config()

def ref():
    """Get the reference sequence"""
    return sml_config['bio.ngs.settings']['db']['ref']
