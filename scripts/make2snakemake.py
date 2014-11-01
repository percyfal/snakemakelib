#!/usr/bin/env python
# File: make2snakemake.py
# Created: Sat Nov  1 13:09:47 2014
# Copyright (C) 2014 by Per Unneberg
#
# Author: Per Unneberg
#

"""
Convert biomake rules to snakemake
"""
import sys
import re
import os

rules = []
config = {}
targets = []
states = [None, "rule", "config"]

header = """
# -*- snakemake -*-
import os
from snakemakelib.utils import update_config, sml_rules_path

# Start by including the general snakefile
include: os.path.join(sml_rules_path(), 'base_settings.rules')
"""

if __name__ == "__main__":
    rule = {}
    if len(sys.argv) != 2:
        print("Usage: ", sys.argv[0], " makefile")
    else:
        state = None
        with open(sys.argv[1], 'r') as fp:
            for l in fp.readlines():
                line = l.rstrip()
                if re.search("^[^ ]+:", line) and not line.startswith("print") and not line.startswith("#"):
                    state = "rule"
                elif line.startswith("ifndef"):
                    state = "config"
                elif line.startswith("endif"):
                    state = None
                elif line.startswith("%"):
                    state = "rule"
                elif line.startswith("\t"):
                    pass
                elif not line:
                    if state == "rule":
                        rules.append(rule)
                    state = None
                elif line.startswith("#"):
                    state = None
                elif re.search("([A-Z_a-z0-9]+)=(.*)", line):
                    if state is None:
                        state = "vardef"
                else:
                    pass
                # Add to config and rule
                if state == "config":
                    if not line.startswith("ifndef") and not line.startswith("ifneq") and not line.startswith("ifeq"):
                        m = re.match("([^=]+)=(.*)", line)
                        (k,v) = m.groups()
                        config[k] = v
                elif state == "rule":
                    print (line)
                    if re.search("^%", line):
                        (output, inputs)  = [x.replace("%", "{prefix}") for x in line.split(":")]
                        rule = {'input': inputs,
                                'output' : output,
                                'shell' : "",
                                'name' : ""}
                    elif re.search("(^[^ %]+):", line):
                        m = re.search("(^[^ %]+):(.*)", line)
                        rule = {'input': "",
                                'output' : "",
                                'shell' : "",
                                'name' : m.group(1)}
                        if m.group(2):
                            rule["input"] = m.group(2)
                    else:
                        rule["shell"] += "{shell}".format(shell=line.lstrip())
                elif state=="vardef":
                    targets.append(line)

        
        print(header)
        # Format config
        print ("config_default = {leftbrace} \n\t\"{section}\" : {leftbrace}".format(leftbrace="{", section=os.path.basename(sys.argv[1].replace("Makefile.", ""))))
        for (section,opt) in config.items():
            print ("\t\t\"{section}\" : \"{opt}\",".format(section=section, opt=opt))
        print ("\t\t{rightbrace},\n\t{rightbrace},\n{rightbrace}".format(rightbrace="}"))

        print ("\nconfig = update_config(config, config_default)\n")

        # print targets
        for tgt in targets:
            print (tgt)
        # Format rules
        i = 0
        for rule in rules:
            i = i + 1
            if rule["name"]:
                print ("rule {name}:".format(name=rule["name"]))
            else:
                print ("rule rule_{cnt}:".format(cnt=i))
            for k in ["input", "output", "shell"]:
                if rule[k]:
                    print ("\t{k}: \"{v}\"".format(k=k, v=rule[k]))


