#!/bin/bash
set -e  # exit on first error

#abinit --version
#abinit --build
abicheck.py --with-flow

nosetests -v --with-coverage --cover-package=pseudo_dojo --logging-level=INFO

# This is to run the integration tests (slow)
# pytest -v --cov=abipy --doctest-modules --ignore=./docs/ abipy/integration_tests
# pytest abipy/integration_tests --ignore=./docs/

# Generate documentation
#if [[ "${ABIPY_SPHINX}" == "yes" ]]; then
#    pip install -r ./docs/requirements.txt
#    ./docs/install_reqs.sh;
#    cd ./docs && make && cd ..
#fi
