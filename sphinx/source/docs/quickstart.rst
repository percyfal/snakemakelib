Quickstart
==========

Introduction
------------

The purpose of snakemakelib is to build a library of rules that can be
reused without actually writing them anew. The motivation is that only
parameters, e.g. program options, inputs and outputs, of a rule change
from time to time, but the rule execution is identical. 

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

2. **Predefined workflows**. The workflows do require a full installation
   of snakemakelib as they depend on functionality therein. See
   :doc:`workflows` for more information.

3. **Simple regexp-based sample organization setup**. Input files and
   samples can be organized in various ways. snakemakelib provides
   regexps for some common sample organization setups, but it is easy
   to add custom configurations. See :doc:`sampleorganization` for
   more information.

4. **Automated generation of alignment indices and annotation file
   names**. Provided a reference is given, snakemakelib offers the
   possibility to automatically generate names of annotation and index
   file names, building the index files if necessary. See
   :doc:`databases` for more information.

   

Quick installation
------------------


Getting started
---------------
