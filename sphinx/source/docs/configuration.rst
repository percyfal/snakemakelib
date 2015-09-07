Configuration
=============

The default configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

For any application, the main rule file contains a *default
configuration* that ensures that all rules have sensible defaults set,
regardless whether the user decides to modify them or not. In
principle, a rule file has two parts:

1. default configuration
2. rules

The basic configuration structure looks like

.. code-block:: text

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

.. code-block:: python

   config_default = {
       'bio.ngs.align.bwa' : {
           'mem' : {
               'options': "",
           },
       },
   }

		
The namespace is ``bio.ngs.align.bwa``, reflecting the fact that the
rules file is located in the folder ``rules/bio/ngs/align`` and is
named ``bwa.rules``. ``config`` is the global snakemake configuration
object. Incidentally, this example shows another key idea of the
configuration, namely that some options inherit from rules files
higher up in the file hierarchy. The rules file
``rules/bio/ngs/settings.rules`` contains a generic configuration that
is common to all ngs rules. This implementation makes it possible to
override settings for specific programs, like for instance the
``threads`` parameter above.

User-defined configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

A user can modify the configuration by defining a :class:`dict` object
and passing it as an argument to
:meth:`~snakemakelib.config.update_config`. This is done in the
Snakefile that uses ``include`` statements to include rules files, and
must be done **before** any ``include`` statement. The reason is that
when a rules file is included, the default configuration values are
compared to the existing ``config``. If the user has defined custom
configurations, these will take precedence over the default values. If
no custom configuration exists, the default values are applied.

As an example, imagine we want to change the number of threads from 8
to 1 for ``bwa mem`` in the example Snakefile above. The modified
Snakefile would then look as follows:

.. code-block:: python

    #-*- snakemake -*-

    # Import config-related stuff
    from snakemakelib.config import update_config
    my_config = {
        'bio.ngs.align.bwa' : {
            'mem' : {
                 'threads' : 1,
            },
        },
    }

    # Update configuration
    update_config(config, my_config)

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
``~/.smlconf.yaml`` if it exists. Typically, here the user would set
variables that are common to all applications, such as paths to
reference databases and commonly used executables. Moreover, if a file
smlconf.yaml is present in the Snakefile working directory, it is loaded
by default. Finally, the user can manually load a file via the
:meth:`~snakemakelib.config.load_sml_config` function.

Environment configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Binary executables that are in the user's PATH environment variable
should be picked up dy default. For some applications, such as java
programs, the search path can be set via an environment variable; for
instance, this is the case for GATK (environment variable GATK\_HOME)
and picard (PICARD\_HOME). Finally, explicit paths can be set in the
configuration file:

.. code-block:: yaml

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

.. code-block:: shell

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

.. code-block:: shell

    biodata/genomes/ORGANISM/BUILD/
        seq/
        star/

In a future version of snakemakelib it will be possible to set the
``bio.ngs.settings.db.build`` variable (e.g. to ``hg19``), provided that
a cloudbiolinux installation is present.
