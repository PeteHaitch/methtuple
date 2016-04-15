#!/usr/bin/env bash
# Based on https://github.com/pysam-developers/pysam/blob/master/run_tests_travis.sh
pushd .

WORKDIR=`pwd`

#Install miniconda python
if [ $TRAVIS_OS_NAME == "osx" ]; then
	curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	bash Miniconda3-latest-MacOSX-x86_64.sh -b
else
	curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b
fi

# Create a new conda environment with the target python version
~/miniconda3/bin/conda install conda-build -y
~/miniconda3/bin/conda create -q -y --name testenv python=$CONDA_PY pysam

# Add new conda environment to PATH
export PATH=~/miniconda3/envs/testenv/bin/:$PATH

# Hack to force linking to anaconda libraries rather than system libraries
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/envs/testenv/lib/
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/miniconda3/envs/testenv/lib/

# Need to make C compiler and linker use the anaconda includes and libraries:
export PREFIX=~/miniconda3/
export CFLAGS="-I${PREFIX}/include -L${PREFIX}/lib"
export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

# create a new folder to store external tools
mkdir -p $WORKDIR/external-tools

popd

# Try building conda recipe first
~/miniconda3/bin/conda-build ci/conda-recipe/ --python=$CONDA_PY

# install code from the repository
python setup.py install
python setup.py test
# TODO: coveralls