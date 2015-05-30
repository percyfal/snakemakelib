The purpose of snakemakelib is to build a library of rules that can be
reused without actually writing them anew. The motivation is that only
parameters, e.g. program options, inputs and outputs, of a rule change
from time to time, but the rule execution is identical. Therefore, my
aim is to provide a very simplistic configuration interface in which the
rule parameters can be modified with simple strings.

Implementation
~~~~~~~~~~~~~~

The backbone configuration object is the class *BaseConfig*, which is a
modified *dict* class. It does simple type checking and also overrides
the *\_\_getitem\_\_*, *\_\_setitem\_\_*, and *update* methods, ensuring
that all entries are *BaseConfig* objects.

An auxiliary function *update\_sml\_config* populates the backend
configuration whenever a include statement loads a rules file. Another
feature is that entries can point to regular python functions.

The default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Each rule file consists of rules and an accompanying *default
configuration*. The latter ensures that all rules have sensible defaults
set, regardless whether the user decides to modify them or not. In
principle, a rule file has two parts:

1. default configuration
2. rules

The basic configuration structure looks like

::

    namespace
        section/parameter
            parameter

The *namespace* is an identifier for the rules file, and should be named
``path.to.rules``, where ``path`` and ``to`` are directory names
relative to the rules root path. The *section/parameter* is either a
parameter related to the program, or a subprogram which in turn can have
*parameters* assigned to it.

As an example, consider the default configuration for
``snakemakelib.bio.ngs.bwa``:

::

    config_default = { 
        'bio.ngs.align.bwa' : {
            'cmd' : "bwa",
            'index' : index,
            'threads' : sml_config['bio.ngs.settings']['threads'],
            'options' : "-M",
            'mem' :{
                'options' : "",
            },
        },
    }

The namespace is ``bio.ngs.align.bwa``, reflecting the fact that the
rules file is located in the folder ``rules/bio/ngs/align`` and is named
``bwa.rules``. ``sml_config`` is a reference to the backend global
*BaseConfig* object that stores all loaded rule configurations.
Incidentally, this example shows another key idea of the configuration,
namely that some options inherit from rules files higher up in the file
hierarchy. The rules file ``rules/bio/ngs/settings.rules`` contains a
generic configuration that is common to all ngs rules. This
implementation makes it possible to override settings for specific
programs, like for instance the ``threads`` parameter above.

Viewing the default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``utils.rules`` defines a rule ``conf`` that can be used to view the
current configuration of included files:

::

    snakemake conf

The output is section according to *namespace*, i.e. the rules file.

Accessing snakemake config
~~~~~~~~~~~~~~~~~~~~~~~~~~

Snakemake defines its own global configuration variable *config* that
can be accessed via the command line. In
``rules/bio/ngs/settings.rules``, three Snakemake *config* options have
been added that are useful in the context of ngs:

::

    # Add configuration variable to snakemake global config object
    config['samples'] = config.get("samples", [])
    config['regions'] = config.get("regions", []) 
    config['runs'] = config.get("runs", [])

For instance, a workflow can be run on a single sample by issuing

::

    snakemake rule --config samples=["SAMPLE1","SAMPLE2"]

User-defined configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

A user can modify the configuration by defining a *BaseConfig* object
and updating the ``sml_config`` object mentioned in the previous
section. This is done in the Snakefile that uses ``include`` statements
to include rules files, and must be done **before** any ``include``
statement. The reason is that when a rules file is included, the default
configuration values are compared to the existing ``sml_config``. If the
user has defined custom configurations, these will take precedence over
the default values. If no custom configuration exists, the default
values are applied.

As an example, imagine we want to change the options to ``-k 10 -w 50``
for ``bwa mem`` in the example Snakefile above. The modified Snakefile
would then look as follows:

::

    #-*- snakemake -*-
    # Snakefile example
    # Add path to snakemakelib, unless installed in a virtualenv or similar
    sys.path.append('/path/to/snakemakelib')

    # Import config-related stuff
    from snakemakelib.config import update_sml_config
    my_config = {
        'bio.ngs.align.bwa' : {
            'mem' : {
                'options' : "-k 10 -w 50",
            },
        },
    }

    # Initialize configuration
    update_sml_config(my_config)

    # Include settings and utilities
    include: "/path/to/snakemake/rules/settings.rules"
    include: "/path/to/snakemake/rules/utils.rules"
    # Include rules for bwa
    include: "/path/to/snakemake/rules/bio/ngs/align/bwa.rules"

Currently, it is necessary to use exactly the same structure as that of
the default configuration for the relevant sections.

Loading configuration from a yaml file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configuration settings can be loaded from external configuration files
in ``yaml`` format. By default, ``snakemakelib`` will load a file
~/.smlconf.yaml if it exists. Typically, here the user would set
variables that are common to all applications, such as paths to
reference databases and commonly used executables. Moreover, if a file
smlconf.yaml is present in the Snakefile working directory, it is loaded
by default. Finally, the user can manually load a file via the
``load_sml_config`` function.

Environment configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Binary executables that are in the user's PATH environment variable
should be picked up dy default. For some applications, such as java
programs, the search path can be set via an environment variable; for
instance, this is the case for GATK (environment variable GATK\_HOME)
and picard (PICARD\_HOME). Finally, explicit paths can be set in the
configuration file:

::

    bio.ngs.align.bwa:
      cmd: /path/to/bwa

The implementation is slightly inconsistent at present. Check the
relevant rules file for what parameter to set.

Configuring reference data
~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, snakemakelib assumes a directory structure for reference
data that is based on the directory layout defined in
`cloudbiolinux <http://cloudbiolinux.org/>`__. Basically, the structure
looks as follows:

::

    biodata/genomes/ORGANISM/BUILD/
        bowtie/
        bowtie2/
        bwa/
        mosaik/
        rnaseq@
        seq/

The ``seq`` subdirectory holds the fasta references, whereas
application-specific indices (e.g. bowtie) are found in a directory with
that name. If you specify the location of the reference sequence in
configuration section ``bio.ngs.settings.db.ref``, snakemakelib will
automatically generate the paths to application-specific indices as
outlined above. In fact, it will even add directories and indices not
defined in cloudbiolinux. To my knowledge, star indices are not
available through cloudbiolinux. Running star will generate indices on
the fly in a directory ``star`` located in the parent directory relative
to the reference file defined in ``bio.ngs.settings.db.ref``. In the
above example, we would have

::

    biodata/genomes/ORGANISM/BUILD/
        seq/
        star/

In a future version of snakemakelib it will be possible to set the
``bio.ngs.settings.db.build`` variable (e.g. to ``hg19``), provided that
a cloudbiolinux installation is present.
