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

snakemakelib is basically a port of the rules in `biomake
<https://github.com/percyfal/biomake>`__ to `Snakemake
<https://bitbucket.org/johanneskoester/snakemake/wiki/Home>`__. The
design principles are similar in that my aim is to compile a library
of rules that can be reused and configured via a simple configuration
interface.


.. warning:: Use the rules at your own risk, and make sure you understand
             them before running any commands. I take no responsibility
             if you'd happen to run a ``snakemake clean`` in an
             inappropriate location, removing precious data in the
             process.

Contents
---------

.. toctree::
   :maxdepth: 2

   docs/quickstart
   docs/installation
   docs/configuration
   docs/sampleorganization
   docs/databases
   docs/troubleshooting


