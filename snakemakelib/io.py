# Copyright (C) 2015 by Per Unneberg
import os
import re
import copy
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def set_output(workflow, temp_rules=[], protected_rules=[], temp_filetypes=[], protected_filetypes=[]):
    """Given list of rules, set their corresponding outputs to temporary or protected.

    Args:
      workflow: snakemake workflow
      temp_rules: list of rules for which output is to be set temporary
      temp_filetypes: file extensions to set to temporary output
      protected_rules: list of rules for which output is to be set protected
      protected_filetypes: file extensions to set to protected output

    Returns:
      None
    """
    _set_output(workflow, temp_rules, temp_filetypes)
    _set_output(workflow, protected_rules, protected_filetypes, label="protected")

def _set_output(workflow, rules, filetypes, label="temp"):
    attr = "{label}_output".format(label=label)
    r = re.compile("|".join("{}$".format(x) for x in list(set(filetypes))))
    for rule in rules:
        try:
            output = copy.copy(workflow._rules[rule].output)
            for o in workflow._rules[rule].output:
                if r.search(o) is None:
                    output.remove(o)
                    smllogger.debug("removing output {output} from {label} output in rule {rule}".format(output=o, rule=rule, label=label))
                if label == "protected" and o in workflow._rules[rule].temp_output:
                    workflow._rules[rule].temp_output.remove(o)
                    smllogger.warn("output {output} present also in temporary output! removing in rule {rule}".format(output=o, rule=rule))
            setattr(workflow._rules[rule], attr, set(output))
            if len(list(set(output))) == 0:
                smllogger.warn("no outputs set to {label} for rule '{rule}'; make sure output file extension is included in list {label}_filetypes; see configuration section 'settings.{label}_filetypes' and 'settings.{label}_filetypes_default'".format(rule=rule, label=label))
        except:
            smllogger.warn("failed to set output in rule {rule} to {label}".format(rule=rule, label=label))
