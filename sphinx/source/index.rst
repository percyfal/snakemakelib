.. snakemakelib documentation master file, created by
   sphinx-quickstart on Fri May 29 14:05:01 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

snakemakelib - a library of snakemake rules
===========================================

.. warning:: snakemakelib is still very much WIP and currently
             undergoing heavy development. I'm aiming for a 0.1.0
             release this autumn.

.. _about:


`Snakemake <https://bitbucket.org/johanneskoester/snakemake/wiki/Home>`__
library for various applications, with a focus on bioinformatics and
next-generation sequencing.

The Snakemake rules contain general recipies for commonly used
applications and bioinformatics programs. The use cases reflect the
needs I've had and do by no means have a comprehensive coverage.
Nevertheless, many commands are so commonly used that the recipes may
be of general interest.

.. warning:: Use the rules at your own risk, and make sure you
             understand them before running any commands. I take no
             responsibility if you'd happen to run a ``snakemake
             clean`` in an inappropriate location, removing precious
             data in the process.


Features
^^^^^^^^

1. **Rule library**. As much as is possible, every rule lives in an
   individual file. This makes it easy to fine-tune what rules to
   include. Moreover, most rules can be included in snakefiles without
   even using functionality in snakemakelib. At the very least, if
   rules need to be tweaked, the rule library serves as a
   cut-and-paste resource of template rules. Rules are located in the
   `rules subfolder
   <https://github.com/percyfal/snakemakelib/tree/master/snakemakelib/rules>`__.

2. **Predefined workflows**. The workflows do require a full
   installation of snakemakelib as they depend on functionality
   therein. See :doc:`docs/workflows` for more information.

3. **Simple regexp-based sample organization setup**. Input files and
   samples can be organized in various ways. snakemakelib provides
   regexps for some common sample organization setups, but it is easy
   to add custom configurations. See :doc:`docs/sampleorganization`
   for more information.

4. **Automated generation of alignment indices and annotation file
   names**. Provided a reference is given, snakemakelib offers the
   possibility to automatically generate names of annotation and index
   file names, building the index files if necessary. See
   :doc:`docs/databases` for more information.


Contents
---------

.. toctree::
   :maxdepth: 2

   docs/quickstart
   docs/installation
   docs/rules
   docs/configuration
   docs/sampleorganization
   docs/databases
   docs/troubleshooting
   docs/release_notes

