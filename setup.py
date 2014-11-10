# Copyright (c) 2014 Per Unneberg
from setuptools import setup, find_packages
import os
import glob

setup(name = "snakemakelib",
      version = "0.1.0",
      author = "Per Unneberg",
      author_email = "per.unneberg@scilifelab.se",
      description = "Snakemake rule library with additional utilities",
      license = "MIT",
      scripts = glob.glob('scripts/*.py'),
      install_requires = [
          "pyyaml",
          "matplotlib>=1.3.1",
          "snakemake>=3.1",
          ## Required for testing
          "nose",
          "sphinx",
      ],
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
