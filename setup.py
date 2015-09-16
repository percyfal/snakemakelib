# Copyright (c) 2014 Per Unneberg
# Modelled on bokeh setup script
# --------------------------------------------------
# Imports
# --------------------------------------------------

from __future__ import print_function

# stdlib
import os
from setuptools import setup
from os.path import realpath, dirname, relpath, join

# Extensions
import versioneer

# --------------------------------------------------
# globals and constants
# --------------------------------------------------

ROOT = dirname(realpath(__file__))

# --------------------------------------------------
# classes and functions
# --------------------------------------------------

package_data = []


def package_path(path, filters=()):
    if not os.path.exists(path):
        raise RuntimeError("packaging non-existent path: %s" % path)
    elif os.path.isfile(path):
        package_data.append(relpath(path, 'snakemakelib'))
    else:
        for path, dirs, files in os.walk(path):
            path = relpath(path, 'snakemakelib')
            for f in files:
                if not filters or f.endswith(filters):
                    package_data.append(join(path, f))

rule_suffixes = ('.rules', '.rule')
workflow_suffixes = ('.workflow')
data_suffixes = ('.fastq.gz', 'fa', '.csv', '.vcf', '.interval_list',
                 '.align_metrics', '.hs_metrics', '.insert_metrics',
                 '.dup_metrics', '.gtf')
                    
package_path(join(ROOT, 'snakemakelib', '_templates'))
package_path(join(ROOT, 'snakemakelib', 'data'))
package_path(join(ROOT, 'snakemakelib', 'rules'), rule_suffixes)
package_path(join(ROOT, 'snakemakelib', 'workflows'), workflow_suffixes)
scripts = []

REQUIRES = [
    'biopython>=1.64',
    'pyyaml>=3.11',
    'snakemake>=3.4.2',
    'texttable>=0.8.2',
    'sphinx>=1.3',
    'pandas>=0.16.0',
    'mock>=1.0.1',
    'pysam>=0.8.3',
    #'bokeh>=0.9.1',
    #'bokehutils==0.1.3',
    'pytest',
    'pytest-cov>=1.8.1',
]

try:
    # Hack for readthedocs
    if not 'readthedocs' in os.path.dirname(os.path.realpath(__file__)):
        REQUIRES.append('matplotlib>=1.4.0')
    else:
        print("readthedocs in path name; assuming we're building docs @readthedocs")
        REQUIRES.append('sphinx-bootstrap-theme')
except:
    pass
    

# https://pythonhosted.org/setuptools/setuptools.html
SETUP_REQUIRES = [
    #'bokehutils==0.1.4',
]    

# Adding github to setup:
# http://mike.zwobble.org/2013/05/adding-git-or-hg-or-svn-dependencies-in-setup-py/
DEPENDENCY_LINKS = [
    #'https://github.com/percyfal/bokehutils/tarball/master@0.1.4#egg=bokehutils'
]

# Integrating pytest with setuptools: see
# https://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
from distutils.core import setup, Command
# you can also import from setuptools

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

_version = versioneer.get_version()
_cmdclass = versioneer.get_cmdclass()
_cmdclass.update({'test': PyTest})
setup(
    name="snakemakelib",
    version=_version,
    cmdclass=_cmdclass,
    author="Per Unneberg",
    author_email="per.unneberg@scilifelab.se",
    description="Snakemake rule library with additional utilities",
    license="MIT",
    url="http://github.com/percyfal/snakemakelib",
    scripts=scripts,
    packages=[
        'snakemakelib',
        'snakemakelib.report',
        'snakemakelib.report.tests',
        'snakemakelib.bio',
        'snakemakelib.bio.ngs',
        'snakemakelib.bio.ngs.align',
        'snakemakelib.bio.ngs.db',
        'snakemakelib.bio.ngs.methylseq',
        'snakemakelib.bio.ngs.qc',
        'snakemakelib.bio.ngs.rnaseq',
        'snakemakelib.bio.ngs.rnaseq.tests',
        'snakemakelib.bio.ngs.tests',
        'snakemakelib.bio.ngs.tools',
        'snakemakelib.bio.ngs.tools.tests',
        'snakemakelib.tests',
        'snakemakelib.workflow',
        'snakemakelib.workflows.tests',
    ],
    package_data={'snakemakelib': package_data},
    setup_requires=SETUP_REQUIRES,
    install_requires=REQUIRES,
    dependency_links=DEPENDENCY_LINKS,
)
