#!/bin/bash

# Setup Environment


workdir=NAME
configfile=NAME.ini
trigtime=1197495364
source /home/maria.okounkova/opt/lscsoft/bayeswave/etc/bayeswave-user-env.sh

bayeswave_pipe --workdir ${workdir} \
    --trigger-time  ${trigtime} \
    --skip-datafind \
    ${configfile} --condor-submit

