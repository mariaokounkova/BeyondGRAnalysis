#!/bin/bash

# Setup Environment


workdir=NAME
configfile=NAME.ini
trigtime=4

bayeswave_pipe --workdir ${workdir} \
    --trigger-time  ${trigtime} \
    --skip-datafind \
    ${configfile} --condor-submit

