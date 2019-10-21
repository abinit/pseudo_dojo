#!/bin/bash
set -e  # exit on first error

abicheck.py --with-flow

#nosetests -v --with-coverage --cover-package=pseudo_dojo --logging-level=INFO

# This is to run the integration tests (slow)
pytest -v --cov=pseudo_dojo --doctest-modules --ignore=./docs/ --ignore=pseudo_dojo/integration_tests

# Generate documentation
#if [[ "${ABIPY_SPHINX}" == "yes" ]]; then
#    pip install -r ./docs/requirements.txt
#    ./docs/install_reqs.sh;
#    cd ./docs && make && cd ..
#fi
