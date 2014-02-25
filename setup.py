import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, 'README.md')).read()
CHANGES = '' # open(os.path.join(here, 'CHANGES.md')).read()

requires = ['pysam']

setup(
  name='comethylation',
  version='0.99.5',
  description='comethylation',
  long_description=README + '\n\n' +  CHANGES,
  classifiers=[
    "Programming Language :: Python",
      "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
  author='Peter Hickey',
  author_email='peter.hickey@gmail.com',
  url='https://github.com/PeteHaitch/Comethylation',
  keywords='bisulfite sequencing methylation bismark',
  packages=find_packages(),
  include_package_data=True,
  zip_safe=False,
  install_requires=requires,
  tests_require=requires,
  test_suite="comethylation.tests",
  scripts = [
    'comethylation/scripts/comethylation',
    'comethylation/scripts/bismarkify',
    'comethylation/scripts/correct_Bismark_PE_SAM'
  ]
)
