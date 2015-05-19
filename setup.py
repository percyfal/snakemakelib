# Copyright (c) 2014 Per Unneberg
from setuptools import setup, find_packages
import os
import glob
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = 'snakemakelib/_version.py'
versioneer.versionfile_build = 'snakemakelib/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = 'snakemakelib-' # dirname like 'myproject-1.2.0'

try:
    with open("requirements.txt", "r") as fh:
        install_requires = [x.strip() for x in fh.readlines()]
except IOError:
    install_requires = []
    
setup(
    name = "snakemakelib",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author = "Per Unneberg",
    author_email = "per.unneberg@scilifelab.se",
    description = "Snakemake rule library with additional utilities",
    license = "MIT",
    url = "http://github.com/percyfal/snakemakelib",
    scripts = glob.glob('scripts/*.py'),
    install_requires = install_requires,
    test_suite = 'nose.collector',
    packages=find_packages(exclude=['ez_setup', 'test*']),
    namespace_packages = [
        'snakemakelib',
    ],
    package_data = {
        'snakemakelib':[
            'data/*',
            'rules/*',
            'examples/*',
            'templates/*',
    ]},
)
