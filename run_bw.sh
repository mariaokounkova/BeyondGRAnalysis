#!/bin/bash

# Setup Environment

workdir=Run
configfile=dCS.ini
trigtime=1126259462.0

bayeswave_pipe --workdir ${workdir} \
    --trigger-time  ${trigtime} \
    --skip-datafind \
    ${configfile} --condor-submit

