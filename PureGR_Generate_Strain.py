#!/usr/bin/python

import sys
sys.path.insert(1, 'catalog_tools')
from convert_sxs_to_lvc import *
import argparse
import scipy
import h5py
import numpy as np
import os

def GetModesFromString(modes):
    """ Method to get output modes from a given string, used to specify
        the modes to the LVC file generation """
    if modes == 'all':
        modes = [[l,m] for l in range(2,9) for m in range(-l,l+1)]
    elif modes == '22only':
        modes = [[2, 2], [2, -2]]
    return modes

def GenerateStrainFiles(only22 = False):

    """ Converts standard rhOverM_Asymptotic_GeometricUnits.h5 file
        into LVC format """


    ## Vars needed for the sxs to lvc conversion script. 
    ## values like resolution and sxs_metadata and sxs_resolutions
    ## don't really matter, but we need to specify them
    sxs_data = "Waveforms" ## input directory
    resolution = 0
    sxs_metadata = "catalog_tools/Metadata/sxs_catalog.json"
    sxs_resolutions = "catalog_tools/Metadata/sxs_catalog_resolutions.json"
    modes = GetModesFromString("22only" if only22 else "all") ## modes to output
    out_path = "Waveforms"
    in_file = "rhOverM_Asymptotic_GeometricUnits.h5"
    out_file = "LVC_Format.h5"

    ## convert the simulation into sxs format
    convert_simulation(sxs_data, resolution, sxs_metadata, sxs_resolutions, modes, out_path, \
                  in_name = in_file, out_name = out_file)
    print("Output LVC format waveform to", out_file)
    
def main():
	p = argparse.ArgumentParser(description="Convert sxs format waveform to lvc format")
	p.add_argument('--only22', help='Only output the 22 mode', \
		dest='only22', action='store_true')
	p.set_defaults(only22=False)
	args = p.parse_args()

	GenerateStrainFiles(args.only22)

if __name__ == "__main__":
  main()
