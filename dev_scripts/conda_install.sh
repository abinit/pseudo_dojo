#!/bin/bash
set -e  # exit on first error

echo "Installing AbiPy dependencies with conda."
echo "Adding conda-forge and abinit to channels"
echo "Working in CONDA_PREFIX: ${CONDA_PREFIX} ..."
conda config --add channels conda-forge

echo "Installing requirements listed requirements.txt and requirements-optional.txt ..."
# https://github.com/ContinuumIO/anaconda-issues/issues/542
conda install -y -c anaconda setuptools
conda install -y --file ./requirements.txt
conda install -y --file ./requirements-optional.txt

#echo "Installing abinit from abinit channel ..."
#conda install -y -c abinit abinit=${ABINIT_VERSION}
#abinit --version
#abinit --build

echo "Installation completed"
