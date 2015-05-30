About
=====

`Snakemake <https://bitbucket.org/johanneskoester/snakemake/wiki/Home>`__
library for various applications, with a focus on bioinformatics and
next-generation sequencing.

The Snakemake rules contain general recipies for commonly used
applications and bioinformatics programs. The use cases reflect the
needs I've had and do by no means have a comprehensive coverage.
Nevertheless, many commands are so commonly used that the recipes may be
of general interest.

Comparison with `biomake <https://github.com/percyfal/biomake>`__\ 
-------------------------------------------------------------------

snakemakelib is basically a port of the rules in
`biomake <https://github.com/percyfal/biomake>`__ to
`Snakemake <https://bitbucket.org/johanneskoester/snakemake/wiki/Home>`__.
The design principles are similar in that my aim is to compile a library
of rules that can be reused and configured via a simple configuration
interface.

Disclaimer
==========

Use the rules at your own risk, and make sure you understand them before
running any commands. I take no responsibility if you'd happen to run a
``snakemake clean`` in an inappropriate location, removing precious data
in the process. You have been warned!

Installation
============

Requirements
------------

-  Snakemake version >= 3.2 that supports the global ``config``
   variable.
-  Python version >= 3.3. This may require setting up a virtual
   environment for python3+.

Obviously you also need working installations of the programs whose
rules you wish to make use of.

Installation instructions
-------------------------

Clone the repository https://github.com/percyfal/snakemakelib to an
appropriate location:

.. code-block:: shell

    git clone https://github.com/percyfal/snakemakelib /path/to/snakemakelib

Then either

-  add the path to your Snakefile, or
-  run ``python setup.py install`` (``python setup.py develop`` for
   development installation)

Running the tests
-----------------

Issue

.. code-block:: shell

    nosetests -v -s

to run a suite of smaller tests. To run a working pipeline, issue

.. code-block:: shell

    nosetests --attr=slow

Quickstart
==========

Probably the best way to get started is to run a real example. See the
`repaper <https://github.com/percyfal/repaper>`__ repository for
workflows that should work out of the box.
