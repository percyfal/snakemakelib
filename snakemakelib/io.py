# Copyright (C) 2015 by Per Unneberg
import os
import re
import copy
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

DEFAULT_TEMP_FILES = ['.bam', '.gz', '.zip', '.bigWig']

def set_temp_output(workflow, rules=[], temp_files=[]):
    """Given list of rules, set their corresponding outputs to temporary outputs.

    Args:
      workflow: snakemake workflow
      rules: list of rules
      temp_files: file extensions 

    Returns:
      None
    """
    sfx = DEFAULT_TEMP_FILES + temp_files
    r = re.compile("|".join("{}$".format(x) for x in sfx))
    for rule in rules:
        smllogger.info("setting output in rule {rule} to temporary".format(rule=rule))
        try:
            output = copy.copy(workflow._rules[rule].output)
            for o in workflow._rules[rule].output:
                if r.search(o) is None:
                    output.remove(o)
            workflow._rules[rule].temp_output = set(output)
        except:
            smllogger.warn("failed to set output in rule {rule} to temporary".format(rule=rule))
