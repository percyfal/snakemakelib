# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, sml_rules_path, BaseConfig, get_sml_config

include: os.path.join(sml_rules_path(), 'settings.rules')
include: os.path.join(sml_rules_path(), 'utils.rules')


include: os.path.join(sml_rules_path(), 'bio/ngs/tools', 'gatk.rules')
include: os.path.join(sml_rules_path(), 'bio/ngs/tools', 'picard.rules')


include: os.path.join(sml_rules_path(), 'bio/ngs/align', 'bwa.rules')
