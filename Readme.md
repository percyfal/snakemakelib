## snakemakelib ##

<table>
<tr>
<td>master</td>	
<td><img src="https://travis-ci.org/percyfal/snakemakelib.svg?branch=master" alt="master" /></td>
<td><img src="https://coveralls.io/repos/percyfal/snakemakelib/badge.svg?branch=master" alt="coverage master" /></td>
</tr>
<tr>
<td>develop</td>
<td><img src="https://travis-ci.org/percyfal/snakemakelib.svg?branch=develop" alt="develop" /></td>
<td><img src="https://coveralls.io/repos/percyfal/snakemakelib/badge.svg?branch=develop" alt="coverage develop" /></td>
</tr>
</table>

## About ##

Snakemake library for all kinds of tasks, with a focus on
bioinformatics and next-generation sequencing. See
[the documentation](http://snakemakelib.readthedocs.org/en/latest/index.html)
for more information.

## Installation ##

snakemakelib is still not on PyPI but you can still install from
github using pip in several ways. I recommend using the
[user installation scheme](https://docs.python.org/3.4/install/index.html#inst-alt-install-user)
by making use of the *--user* flag. You need to set the environment
variable
[PYTHONUSERBASE](https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONUSERBASE)
to e.g. *$HOME/lib/python*. Then, with this setting, set
[PYTHONPATH](https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONPATH)
to *$HOME/lib/python/lib/python3.4/site-packages* to access the
package. The examples below show the user installation scheme.

Alternatively, install the package into a
[virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/).

### Installing directly from github ###

You can install as an
[editable install](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs)
by invoking

	pip3 install git+https://github.com/percyfal/snakemakelib.git@master#egg=snakemakelib --user

### Cloning and installing ###

The previous command currently fails to install
[bokehutils](https://github.com/percyfal/bokehutils) (see
[issue 16](https://github.com/percyfal/snakemakelib/issues/16)). As a
workaround, you can manually clone the repo

	git clone git@github.com:percyfal/snakemakelib.git

and invoke (standing in the cloned repo)

	pip3 install -r requirements.txt --user

## Quickstart ##

Probably the best way to get started is to run a real example. See the
[repaper](<https://github.com/percyfal/repaper>) repository for
workflows that should work out of the box.
