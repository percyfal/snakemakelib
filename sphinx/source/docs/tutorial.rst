Tutorial
========

The intended usage is that the user first creates a Snakefile for use
with a particular dataset/problem. Thereafter, include statements are
used to include rules of interest.

The tutorial Snakefiles can be found in the ``examples`` directory.

A simple example
----------------

Let's start by writing a Snakefile that includes the rules for the
aligner ``bwa``.

.. code-block:: python

    #-*- snakemake -*-
    # Snakefile example

    from snakemakelib.config import sml_rules_path

    # Include settings and utilities
    include: os.path.join(sml_rules_path(), "utils.rules")
    include: os.path.join(sml_rules_path(), "bio/ngs/settings.rules")
    # Include rules for bwa
    include: os.path.join(sml_rules_path(), "bio/ngs/align/bwa.rules")

Here, we make use of a function ``sml_rules_path`` which dynamically
calculates the absolute installation path of the ``snakemakelib/rules``
directory. Snakemake includes options to view tasks:

.. code-block:: shell

    $ snakemake -l

    conf
        Print global configuration settings for loaded rules
    conf_sections
        Print configuration sections. TODO: insert docstring from each rule file without running an import statement
    rule_ll
        Print rules sorted by rule file, including input/output information
    rule_l
        Print rules sorted by rule file, excluding input/output information
    samtools_sam2bam
        Convert sam file to bam
    samtools_gene_region_from_bam
        Extract gene region from bam file
    samtools_index
        Run samtools index
    bwa_mem
        Run bwa mem

The ``bwa_mem`` rule comes from ``bwa.rules``. Note that also rules for
samtools are present; these have been included via an include statement
in ``bwa.rules``. In addition, the rules file ``utils.rules`` included
above defines convenience rules for viewing more detailed rule
information. For instance, the rule ``rule_ll`` shows the following:

.. code-block:: shell

    $ snakemake rule_ll

    Provided cores: 1
    Job counts:
        count   jobs
        1   rule_ll
        1
    rule rule_ll:
    utils
    =====
    conf

        Print global configuration settings for loaded rules
    conf_sections

        Print configuration sections. TODO: insert docstring from each rule file without running an import statement
    rule_l

        Print rules sorted by rule file, excluding input/output information
    rule_ll

        Print rules sorted by rule file, including input/output information
    bio.ngs.align.bwa
    =================
    bwa_mem

        Run bwa mem
        input: {prefix}_1.fastq.gz {prefix}_2.fastq.gz
        output: {prefix}.bam
        shell: bwa mem -t {threads}  bwa/ {prefix}_1.fastq.gz {prefix}_2.fastq.gz | samtools view -Sb - > {prefix}.bam

    bio.ngs.tools.samtools
    ======================
    samtools_gene_region_from_bam

        Extract gene region from bam file
        input: {prefix}.bam {prefix}.region_{region}.{sfx}.bed
        output: {prefix}.region_{region}.{sfx}.bam
        shell: samtools view  -b -L {prefix}.region_{region}.{sfx}.bed {prefix}.bam > {prefix}.region_{region}.{sfx}.bam

    samtools_index

        Run samtools index
        input: {prefix}.bam
        output: {prefix}.bai
        shell: samtools  index {prefix}.bam

    samtools_sam2bam

        Convert sam file to bam
        input: {prefix}.sam
        output: {prefix}.bam
        shell: samtools view  -Sb {prefix}.sam > {prefix}.bam

    1 of 1 steps (100%) done

This rule prints the docstring and the definitions of the ``input``,
``output``, and ``shell`` parameters. In the example above, we see that
the output for rule ``bwa_mem`` looks like ``{prefix}.bam``, where
``{prefix}`` is a wildcard that matches a given pattern in the inputs.
To see what they look like in this example, run

.. code-block:: shell

    $ snakemake test.bam

    MissingInputException in line 37 of /path/to/snakemakelib/bio/ngs/align/bwa.rules:
    Missing input files for rule bwa_mem:
    test_1.fastq.gz
    test_2.fastq.gz

As its name implies, Snakemake works like GNU Make in that one seeks to
build a target output, in this case ``test.bam``. Had the files
``test_1.fastq.gz`` and ``test_2.fastq.gz`` been present, the rule
``bwa_mem`` had run the command defined in the ``shell`` section:

::

    bwa mem -t {threads}  bwa/ test_1.fastq.gz test_2.fastq.gz | samtools view -Sb - > test.bam

where ``{threads}`` will be substituted for the number of threads
configured for bwa. Note that the index is ``bwa/`` as no reference has
yet been configured. We'll see in the next example how these can be
configured. For now, we can have a look at what is configured as
follows:

.. code-block:: shell

    $ snakemake conf

    Snakemake global configuration variables
    ========================================
     regions : []
     runs : []
     samples : []
     section : 

    Top-level configuration variables
    =================================

    bio.ngs.settings
    ================
        aligner : bwa
        annotation : {'annot_label': '', 'transcript_annot_gtf': ''}
        center : 
        db : {'dbsnp': '', 'ref': '', 'build_config': None, 'build': ''}
        fastq_suffix : .fastq.gz
        inputdir : .
        java : {'java_tmpdir': '/tmp', 'java_mem': '6g'}
        read1_label : _1
        read1_suffix : .fastq.gz
        read2_label : _2
        read2_suffix : .fastq.gz
        regions : []
        rnaseq : {'quantification': ['rsem']}
        runs : []
        sample_column_map : {}
        sample_organization : sample
        sampleinfo : 
        sampleorg : sample_organization(raw_run_re={}, run_id_re={}, sample_re={})
        samples : []
        sequence_capture : {'bait_regions': '', 'target_regions': ''}
        threads : 8
    bio.ngs.align.bwa
    =================
        cmd : bwa
        index : <function index at 0x7fe523d37158>
        mem : {'options': ''}
        options : -M
        threads : 8
    bio.ngs.tools.samtools
    ======================
        cmd : samtools
        index : {'options': ''}
        options : 
        ref : 
        threads : 8

A slightly more complicated example
-----------------------------------

This example is partly modelled on how I set up the `ATAC-seq protocol
<https://github.com/percyfal/snakemakelib/blob/develop/workflows/bio/ATAC-seq.workflow>`__

Importing relevant libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start by importing relevant libraries:

::

    # -*- snakemake -*-
    import os
    import pysam
    from snakemakelib.io import set_temp_output
    from snakemakelib.config import update_sml_config, get_sml_config, sml_rules_path
    from snakemakelib.bio.ngs.targets import generic_target_generator

``set_temp_output`` is a helper function for dynamically setting the
output of given rules to temporary. ``update_sml_config`` and
``get_sml_config`` are two important functions; they communicate with
the backend configuration, i.e. the global state of the system.
Configuration options are set and retrieved via these functions.
``sml_rules_path`` retrieves the path to the root location of
snakemakelib rules. Finally, ``generic_target_generator`` is a function
that will generate target names based on the input data. The results of
this function is a list of file names that serve as targets for the
snakemake workflow.

Setting the default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The next step is to set default configuration values for the workflow.
Here I've included a subsection of these defaults:

.. code-block:: python

    # Default configuration settings custom-tailored for ATAC-Seq analysis
    atac_config = {
        'workflows.bio.atac_seq' : {
            'aligner' : 'bowtie',
            'peakcallers' : ['zinba', 'dfilter', 'macs2'],
            'trimadaptor' : True,
        },
        'settings' : {
            'temp_rules' : [],
        },
        'bio.ngs.enrichment.macs' : {
            'callpeak' : {
                'options' : '-g dm --nomodel --shiftsize 50 -q 0.01',
            },
        },
    }

To begin with the obvious the configuration object is a dictionary of
dictionaries. The top-level keys correspond to configuration sections,
and for each step we descend in the hierarchy, the more specific the
configuration settings are. The section names reflect the locations of
the corresponding rules files. For instance, the rules for section
``settings`` is found in ``rules/settings.rules`` relative to the
snakemakelib installation directory. Similarly, the rules file for
section ``bio.ngs.enrichment.macs`` is located at
``rules/bio/ngs/enrichment/macs.rules``. Note that these files contain a
more comprehensive list of configuration default. The purpose of setting
them in the workflow file is to ``override`` the default settings so
that they are fine-tuned to this particular analysis.

The configuration section ``workflows.bio.atac_seq`` relates to the
workflow file itself. The intention of the parameters in this section is
to govern the behaviour of the workflow in general. For instance, in
this example, we have settings for choosing aligner and what peakcallers
to use.

Finally, it is extremely important to make these configuration settings
visible to the snakemakelib configuration backend. This is done with the
following code:

::

    update_sml_config(atac_config)

Making use of a configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As noted above, the configuration settings override the default settings
in the rules files. However, by supplying a yaml configuration file it
is possible to override these settings too. In `the Buenrostro repaper
example <https://github.com/percyfal/repaper/tree/master/Buenrostro_2013>`__,
the Snakefile contains the following code:

::

    # Load external configuration files
    config = load_sml_config(config)

The documentation for :meth:`~snakemakelib.config.load_sml_config` explains what it does:

.. autofunction:: snakemakelib.config.load_sml_config

Consequently, one can have several layers of configurations, tailored
for different analyses or computing environments.

A digression on sample organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before we delve into the task of generating target names for snakemake,
a digression on sample organization is in order. In the `Buenrostro
configuration
file <https://github.com/percyfal/repaper/tree/master/Buenrostro_2013/smlconf.yaml>`__,
there configuration section ``bio.ngs.settings`` reads as follows:

.. code-block:: yaml

    bio.ngs.settings:
      sample_organization: sample_run_sra

The key ``sample_organization`` tells snakemakelib how data is and shall
be organized. Typically, sample data is of three types:

1. raw run data, i.e raw data from the sequencing machine
2. processed run data
3. sample-level data in which data from several runs for a sample have
   been merged

Since snakemake works on file names, here, we're interested in the
format of file names. As we'll see, snakemakelib determines these names
based on regular expressions. Setting the key ``sample_organization`` to
``sample_run_sra`` will use a predefined set of regular expressions for
the three levels of sample file names described above. In
``bio.ngs.settings``, there is a dictionary ``sample_organization`` that
defines regular expressions for different naming conventions. For the
``sample_run_sra`` key, we have

.. code-block:: python

    sample_org = namedtuple('sample_organization', 'raw_run_re run_id_re sample_re')
    sample_organization = {
        # Data in sample directory divided in subdirectory for each run, SRA-like
        'sample_run_sra' : {
            'sampleorg' : sample_org(RunRegexp(os.path.join("(?P<SM>[a-zA-Z0-9]+)", "(?P<PU>[a-zA-Z0-9]+)", "(?P=PU)")),
                                     RunRegexp(os.path.join("(?P<SM>[a-zA-Z0-9]+)", "(?P<PU>[a-zA-Z0-9]+)", "(?P=PU)")),
                                     SampleRegexp(os.path.join("(?P<SM>[a-zA-Z0-9]+)", "(?P=SM)"))),
            },

Setting ``sample_organization`` to ``sample_run_sra`` will load a
dictionary with key ``sampleorg`` whose value is a :class:`namedtuple`
consisting of three regular expression classes. The regular
expressions use *symbolic group names* (here ``SM`` and ``PU``) to
access the matches in each parenthesis. SRA run names are on the
format ``SRR######`` (see `Understanding SRA Search Results
<http://www.ncbi.nlm.nih.gov/books/NBK56913/#search.what_do_the_different_sra_accessi>`__),
sample names ``SRS######``, where ``#`` represents a digit. As can be
seen above, snakemakelib assumes the following format for the three
levels (even though there is no specific requirement that they begin
with SRS and SRA in the regular expression):

1. raw_run:  ``SRS###### / SRA###### / SRA######``
2. run_id: ``SRS###### / SRA###### / SRA######``
3. sample: ``SRS###### / SRS######``

The symbolic group names ``SM`` and ``PU`` are not randomly chosen;
the same tags are used here as for read group information in the `SAM
format specification
<https://www.google.se/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0CCEQFjAA&url=https%3A%2F%2Fsamtools.github.io%2Fhts-specs%2FSAMv1.pdf&ei=Yyo6VaqIIYSMsAGBqoG4DA&usg=AFQjCNHFmjxTXKnxYqN0WpIFjZNylwPm0Q&bvm=bv.91427555,d.bGg>`__.
Snakemakelib will use this information to write readgroup information
to bam headers, if requested.

A more complicated example is shown below:

.. code-block:: python

    sample_organization = {
        ...
        # Illumina sequence data as delivered by SciLife
        'Illumina@SciLife' : {
            'sampleorg' : sample_org(RunRegexp(os.path.join("(?P<SM>P[0-9]+_[0-9]+)", "(?P<DT>[0-9]+)_(?P<PU1>[A-Z0-9]+XX)", "(?P<PU2>[0-9])_(?P=DT)_(?P=PU1)_(?P=SM)")),
                                     RunRegexp(os.path.join("(?P<SM>P[0-9]+_[0-9]+)", "(?P<DT>[0-9]+)_(?P<PU1>[A-Z0-9]+XX)", "(?P<PU2>[0-9])_(?P=DT)_(?P=PU1)_(?P=SM)")),
                                     SampleRegexp(os.path.join("(?P<SM>P[0-9]+_[0-9]+)", "(?P=SM)"))),
            },

The `example data <../data/projects>`__ are formatted according to this
convention. For instance, for project J.Doe\_00\_01 the full path to the
fastq file is

::

    P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101_1.fastq.gz

Matching this up with the raw run regular expression above, we get
``SM=P001_101``, ``DT=120924``, ``PU1=AC003CCC``, and ``PU2=1``.
``PU`` stands for Platform Unit and could be 'flowcell-barcode.lane'
for Illumina data. This information is split above, hence the
indexing. When snakemakelib writes the ``PU`` read group tag, it will
concatenate the indexed ``PU`` tags in successive order with an
underscore, so that ``PU=PU1_PU2`` above would be ``AC003CCC_1``.
Hence, as sequencing data is a common starting point, given the file
organization and regular expressions above, snakemakelib will
automatically generate names for downstream processing.

Why the need for raw run versus run identifiers? Well, it turns out that
some sequencing centers deliver sequencing data in one format, whereas
the analysis pipelines transform the names into another format.
Previously, I worked with data where the raw data file names looked as
above, but where a processed run would look like

::

    P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_P001_101_1.fastq.gz

but the processed run data prefix would look like

::

    P001_101/120924_AC003CCCXX/1_120924_AC003CCCXX_1

Having two levels of regular expressions for run data allows for this
flexible treatment of file names.

Using a samplesheet
^^^^^^^^^^^^^^^^^^^

NB: the following has only currently been tested for
``bio.ngs.tools.sratools``

Finally, there is another setting in ``bio.ngs.settings`` in the
`Buenrostro configuration
file <https://github.com/percyfal/repaper/tree/master/Buenrostro_2013/smlconf.yaml>`__:

.. code-block:: yaml

    sample_column_map:
        SampleName: SM
        Run: PU

If a samplesheet is present, snakemakelib will use this file instead of
raw file names. ``sample_column_map`` maps the column names in the
samplesheet file to read group identifiers. You must at least provide a
mapping for the read group identifiers present in the regular
expressions. The column names above are from the SRA project information
file. Incidentally, for expression experiments, the ``SampleName``
actually holds GEO identifiers, so example paths to the three sample
file name levels in the Buenrostro example would look like

::

    raw_run: GSM1155957/SRR891268/SRR891268_1.fastq.gz
    run_id:  GSM1155957/SRR891268/SRR891268.bam
    sample:  GSM1155957/GSM1155957.sort.merge.bam

Generating targets
~~~~~~~~~~~~~~~~~~

Writing a custom rule
~~~~~~~~~~~~~~~~~~~~~

