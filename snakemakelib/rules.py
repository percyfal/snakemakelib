# Copyright (C) 2015 by Per Unneberg
import os
import copy
from snakemake.workflow import Workflow
from snakemake.io import OutputFiles
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def create_rule_from_existing(name, template, workflow, **kw):
    """Create rule from a template.
    
    Create rule from existing rule and add it to the workflow. By
    passing keyword arguments it is also possible to update/modify the
    input, output and/or params.

    Args:
        name (str): name of new rule
        template (str): name of existing template rule
        workflow (:class:`Workflow <snakemake.workflow.Workflow>`): snakemake workflow
        kw (dict): keyword argument for updating input, output and/or params
        
    Returns:
        None
    """
    assert type(workflow) is  Workflow, "workflow is not a Workflow: {}".format(workflow)
    try:
        rule = copy.copy(workflow.get_rule(template))
        rule.name = name
    except:
        smllogger.warn("no such template rule '{}'; make sure you have included the template rule file".format(template))
        raise
    workflow.add_rule(name=rule.name)
    workflow._rules[name] = rule
    if kw.get('output'):
        assert type(kw['output']) is tuple, "output argument must be a tuple of type (tuple, dict)"
        rule._output = OutputFiles()
        workflow._rules[name].set_output(*kw['output'][0], **kw['output'][1])
    if kw.get('input'):
        assert type(kw['input']) is tuple, "input argument must be a tuple of type (tuple, dict)"
        workflow._rules[name].set_input(*kw['input'][0], **kw['input'][1])
    if kw.get('params'):
        assert type(kw['params']) is tuple, "params argument must be a tuple of type (tuple, dict)"
        workflow._rules[name].set_params(*kw['params'][0], **kw['params'][1])
