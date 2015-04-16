# Copyright (c) 2014 Per Unneberg
from setuptools import setup, find_packages
import os
import glob

try:
    with open("requirements.txt", "r") as fh:
        install_requires = [x.strip() for x in fh.readlines()]
except IOError:
    install_requires = []
    
setup(name = "snakemakelib",
      version = "0.1.0",
      author = "Per Unneberg",
      author_email = "per.unneberg@scilifelab.se",
      description = "Snakemake rule library with additional utilities",
      license = "MIT",
      scripts = glob.glob('scripts/*.py'),
      install_requires = install_requires,
      test_suite = 'nose.collector',
      packages=find_packages(exclude=['ez_setup', 'test*']),
      namespace_packages = [
          'snakemakelib',
          'snakemakelib.ext',
      ],
      package_data = {
          'snakemakelib':[
              'rules/*',
          ]},
  )

os.system("git rev-parse --short --verify HEAD > ~/.snakemakelib_version")
