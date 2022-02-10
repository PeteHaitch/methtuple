import os

from setuptools import setup, find_packages

# Convert README.md to README.rst automagically
# https://packaging.python.org/en/latest/guides/making-a-pypi-friendly-readme/
# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
  name='methtuple',
  version='1.7.0',
  description='methtuple',
  long_description=long_description,
  long_description_content_type='text/markdown',
  license='MIT',
  classifiers=[
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    "License :: OSI Approved :: MIT License",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
  author='Peter Hickey',
  author_email='peter.hickey@gmail.com',
  url='https://github.com/PeteHaitch/methtuple',
  keywords='bisulfite sequencing methylation bismark bioinformatics',
  packages=find_packages(),
  include_package_data=True,
  zip_safe=False,
  install_requires=['pysam >= 0.8.4'],
  test_suite="methtuple.tests",
  scripts = [
    'methtuple/scripts/methtuple'
  ]
)
