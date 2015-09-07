Quickstart
==========

Introduction
------------

The purpose of snakemakelib is to build a library of rules that can be
reused without actually writing them anew. The motivation is that only
parameters, e.g. program options, inputs and outputs, of a rule change
from time to time, but the rule execution is identical. 


Quick installation
------------------

Install directly from github with pip using the following command:

.. code-block:: shell
		
   pip3 install -e git+https://github.com/percyfal/snakemakelib.git@master#egg=snakemakelib --user

See :doc:`installation` for more information.   

Getting started
---------------

Most importantly, snakemakelib offers a library of rules that can be
included in a snakefile and configured via an external configuration
file.

Bwa alignment
^^^^^^^^^^^^^^

snakemakelib comes packaged with a small `test data set
<https://github.com/percyfal/snakemakelib/tree/master/data>`_
including input data and reference sequences. Provided you have `bwa
<http://bio-bwa.sourceforge.net/>`_ installed, you should be able to
setup the alignment rules as follows.

First, create a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
with the following content:

.. code:: python

   # -*- snakemake -*-
   import os
   from snakemakelib.config import update_config, SNAKEMAKELIB_PATH

   workdir: os.path.join(SNAKEMAKELIB_PATH, "data/projects/J.Doe_00_01"

   include: os.path.join(SNAKEMAKELIB_PATH, "rules/bio/ngs/align/bwa.rules")

   local_config = {
	  '
   }

Variant calling
^^^^^^^^^^^^^^^^

Running the tests
-----------------
