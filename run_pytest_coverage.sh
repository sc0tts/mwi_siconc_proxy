#!/bin/bash
#
# ./run_install_test_cover.sh
#
# Runs pytest
# Then, if pytest runs without error, then run coverage test

pytest; retval=$? ; if [ "$retval" == "0" ]; then ./tests/covtest ; fi
