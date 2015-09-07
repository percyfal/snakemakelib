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
included in a Snakefile and configured via an external configuration
file. Snakemake works in a manner similar to `GNU Make
<https://www.gnu.org/software/make/>`_ in that rules determine how to
generate output files from input files. Output file names are matched
against input file names, whereby wildcards can be used to write
general rules. This feature has been adopted heavily in snakemakelib.
In fact, most rules are of the form

.. code:: python

   rule somerule:
       input: "{prefix}.inputsuffix"
       output: "{prefix}.outputsuffix"
       run: "command {input} > {output}"

where ``{prefix}`` denotes a wildcard. ``.inputsuffix`` and
``.outputsuffix`` are generally application dependent. For instance,
in the bwa example that follows, the input suffix is generally
``.fastq.gz`` and the output suffix ``.bam``.

Bwa alignment
^^^^^^^^^^^^^^

.. note:: This example requires you have installed `bwa
          <http://bio-bwa.sourceforge.net/>`_ and `samtools
          <http://www.htslib.org/>`_

snakemakelib comes packaged with a small `test data set
<https://github.com/percyfal/snakemakelib/tree/master/data>`_
including input data and reference sequences. 

First, create a `Snakefile
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-writing-snakefiles>`_
with the following content:

.. code:: python

   # -*- snakemake -*-
   import os
   from snakemakelib.config import SNAKEMAKELIB_PATH
   
   config = {
       'bio.ngs.align.bwa' : {
	   'index' : os.path.join(SNAKEMAKELIB_PATH, "data/genomes/Hsapiens/hg19/bwa/chr11.fa"),
       },
   }

   workdir: os.path.join(SNAKEMAKELIB_PATH, "data/projects/J.Doe_00_01")
   include: os.path.join(SNAKEMAKELIB_PATH, "rules/bio/ngs/align/bwa.rules")

Briefly, in this file we:

1. set a **configuration dictionary** that indicates where the bwa
   index files are located
2. set a **working directory**
3. **include** rules for bwa

In this minimal example, we make use of the internal snakemakelib
variable ``SNAKEMAKELIB_PATH`` to locate the rules, but could as well
have included the rules by supplying the full path.

Now, to see which rules are included, you can type:

.. code:: shell

   $ snakemake -l

which should generate the following output:
   
.. code:: shell

   bwa_index
	bwa index a reference
   samtools_index
	Run samtools index
   bwa_mem
	Run bwa mem

So, by including ``bwa.rules``, we have actually defined three
`snakemake rules
<https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-rules>`_.

Now, in the test data set relative to the working directory defined
above there are sequence input files organized under samples and
sequencing runs. For instance, in subfolder
``P001_101/120924_AC003CCCXX`` we have the files

.. code:: shell

   1_120924_AC003CCCXX_P001_101_1.fastq.gz
   1_120924_AC003CCCXX_P001_101_2.fastq.gz

Since the bwa rule has the output suffix ``.bam`` and input suffix
``.fastq.gz``, we can align these input files by issuing [#f1]_

.. code:: shell

   snakemake -F P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101.bam

In addition to performing the alignment, this command will generate
bwa indices on the fly. The flag ``-F`` tells snakemake to rerun the
rules, even if outputs are present.
	

Running the workflow tests
--------------------------

The above example is included in the workflow test suite and can be
run with the following command:

.. code:: shell

   py.test -vsk test_bwa_align

There is also a larger example that runs a variant calling pipeline.
In addition to bwa and samtools, this example requires you have `GATK
<https://www.broadinstitute.org/gatk/>`_ and `picard
<http://broadinstitute.github.io/picard/>`_ installed, and you need to
set the environment variables ``GATK_HOME`` and ``PICARD_HOME`` to
point to the installation directories. With this done, you can run the
variant calling pipeline like so:

.. code:: shell

   py.test -m slow
   

.. rubric:: Footnotes

.. [#f1] Actually, the rule also takes into account the read label(s),
         here ``_1`` and ``_2``. These are omitted in the output name.
