#!/bin/bash

conda deactivate
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/etc/profile.d/conda.sh
conda activate igwn-py38
export LAL_DATA_PATH=/home/maria.okounkova/.local/lib/python3.6/site-packages/gwsurrogate/surrogate_downloads/:$LAL_DATA_PATH
python3 /home/maria.okounkova/BeyondGRAnalysis/BilbyPE/FRAMES_DIR/Frames.py
