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

-  Snakemake version >= 3.3 that supports the global ``config``
   variable.
-  Python version >= 3.3. This may require setting up a virtual
   environment for python3+.

Obviously you also need working installations of the programs whose
rules you wish to make use of.

Installation instructions
-------------------------

snakemakelib is still not on PyPI but you can still install from
github using pip in several ways. I recommend using the `user
installation scheme
<https://docs.python.org/3.4/install/index.html#inst-alt-install-user>`_
by making use of the *--user* flag. You need to set the environment
variable `PYTHONUSERBASE
<https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONUSERBASE>`_
to e.g. *$HOME/lib/python*. Then, with this setting, set `PYTHONPATH
<https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONPATH>`_
to *$HOME/lib/python/lib/python3.4/site-packages* to access the
package. The examples below show the user installation scheme.

Alternatively, install the package into a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.

Installing directly from github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


You can install as an `editable install
<https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs>`_
by invoking

.. code-block:: shell

    pip3 install git+https://github.com/percyfal/snakemakelib.git@master#egg=snakemakelib --user


Cloning and installing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous command currently fails to install `bokehutils
<https://github.com/percyfal/bokehutils>`_ (see `issue 16
<https://github.com/percyfal/snakemakelib/issues/16>`_). As a
workaround, you can manually clone the repo

.. code-block:: shell

    git clone git@github.com:percyfal/snakemakelib.git

and invoke (standing in the cloned repo)

.. code-block:: shell

    pip3 install -r requirements.txt --user .

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
