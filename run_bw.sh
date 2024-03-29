#!/bin/bash

# Setup Environment

workdir=BW_${PWD##*/}
configfile=dCS.ini
trigtime=1126259462.4107006

bayeswave_pipe --workdir ${workdir} \
    --trigger-time  ${trigtime} \
    --skip-datafind \
    ${configfile} --condor-submit

