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

The user install scheme
-----------------------

snakemakelib is still not on PyPI but you can still install from
github using pip in several ways. I recommend using the `user
installation scheme
<https://docs.python.org/3.4/install/index.html#inst-alt-install-user>`__
by making use of the *--user* flag.

By default, user site-packages are placed in *~/.local*. If you choose
to install in a different location, you need to set the environment
variable `PYTHONUSERBASE
<https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONUSERBASE>`__
to e.g. *$HOME/lib/python*. Then, with this setting, set `PYTHONPATH
<https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONPATH>`__
to *$HOME/lib/python/lib/python3.4/site-packages* to access the
package. The examples below show the user installation scheme.

Alternatively, install the package into a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/)>`__.

Dependencies
------------

Installation directly from github unfortunately does not install the
dependency `bokehutils <https://github.com/percyfal/bokehutils>`__
(see `issue 16
<https://github.com/percyfal/snakemakelib/issues/16>`__). Either
install bokehutils

.. code-block:: shell

   pip3 install -e git+https://github.com/percyfal/bokehutils.git@develop#egg=bokehutils --user


Installation from github with pip
---------------------------------

You can install as an `editable install
<https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs>`__
by invoking

.. code-block:: shell
		
   pip3 install -e git+https://github.com/percyfal/snakemakelib.git@master#egg=snakemakelib --user

This will install the latest stable version. To stay on bleeding edge,
you can install the develop branch

.. code-block:: shell
		
   pip3 install -e git+https://github.com/percyfal/snakemakelib.git@develop#egg=snakemakelib --user



Installing from source
-----------------------

Alternatively, you can install from source by cloning the snakemakelib
repo

.. code-block:: shell

   git clone git@github.com:percyfal/snakemakelib.git

and invoke (standing in the cloned repo)

.. code-block:: shell
		
   pip3 install -r requirements.txt --user -e .
