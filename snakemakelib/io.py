# Copyright (C) 2015 by Per Unneberg
import os
import re
import copy
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def set_temp_output(workflow, rules, temp_filetypes):
    """Given list of rules, set their corresponding outputs to temporary outputs.

    Args:
      workflow: snakemake workflow
      rules: list of rules
      temp_filetypes: file extensions to set to temporary output

    Returns:
      None
    """
    r = re.compile("|".join("{}$".format(x) for x in list(set(temp_filetypes))))
    for rule in rules:
        smllogger.info("setting output in rule {rule} to temporary".format(rule=rule))
        try:
            output = copy.copy(workflow._rules[rule].output)
            for o in workflow._rules[rule].output:
                if r.search(o) is None:
                    output.remove(o)
                    smllogger.debug("removing output {output} from temporary output in rule {rule}".format(output=o, rule=rule))
            workflow._rules[rule].temp_output = set(output)
        except:
            smllogger.warn("failed to set output in rule {rule} to temporary".format(rule=rule))
