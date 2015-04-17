# Copyright (C) 2015 by Per Unneberg
import os
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def set_temp_output(workflow, rules=[]):
    """Given list of rules, set their corresponding outputs to temporary outputs"""
    for rule in rules:
        smllogger.info("setting output in rule {rule} to temporary".format(rule=rule))
        try:
            workflow._rules[rule].temp_output = set(workflow._rules[rule].output)
        except:
            smllogger.warn("failed to set output in rule {rule} to temporary".format(rule=rule))
