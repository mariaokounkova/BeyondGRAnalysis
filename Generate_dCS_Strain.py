#!/usr/bin/python

import sys
import argparse
import scipy
import h5py
import numpy as np
import os
sys.path.insert(1, 'catalog_tools')
from convert_sxs_to_lvc import *

def ReadExtrapolatedMode(p, piece, mode, lev, order=2, ell=None):
  """ Given a file of extrapolated modes, read in the (mode)
      at a given order """
  ell = str(ell).replace('.', 'p')
  piece_dict = {"DeltaStrain" : "Lev" + str(lev) + "/DeltaStrain.h5", \
          "BackgroundStrain" : "Lev" + str(lev) + "/BackgroundStrain.h5", \
          "dCSModified" : "/dCS_" + ell + "_Lev" + str(lev) + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5", \
                  "hRWZ" : "Lev" + str(lev) + "/rhOverM_Asymptotic_GeometricUnits.h5"}
  file = p + piece_dict[piece]
  l = mode[0]
  m = mode[1]
  f = h5py.File(file, 'r')
  data = f['Extrapolated_N'+str(order)+'.dir']['Y_l' + str(l) + '_m'  + str(m) + '.dat']
  time, re, im = data[:,0], data[:,1], data[:,2]
  result = re + 1j*im
  return time, result

def ComputedCSModifiedStrain(p, mode, l, lev):
  """ Given a value of the dCS coupling constant l, a path 
      p to the extrapolated hPsi4 compute the modified gravitational wave strain """

  ## Read in the background
  time, strain = ReadExtrapolatedMode(p, "BackgroundStrain", mode, lev)
  delta_time, delta_strain = ReadExtrapolatedMode(p, "DeltaStrain", mode, lev)

  ## Now add the strain and delta strain together
  ## with the correct value of l
  total = strain + l**4 * delta_strain

  return time, total


def OutputdCSModifiedStrain(p, ell, only22, dropm0, lev):
    """ Generate an h5 file with the modified strain for a 
        given value of ell """

    ## For naming the file, replace . with p because otherwise
    ## the .h5 file can't be read by catalog scripts
    name = str(ell).replace('.', 'p')
    name_lev = str(ell).replace('.', 'p') + '_Lev' + str(lev)
    
    ## Make the output directory
    os.mkdir(p + '/dCS_' + name_lev)
    
    OutFile = p + 'dCS_' + name_lev + \
      '/rhOverM_Asymptotic_GeometricUnits_CoM.h5'
    fOut = h5py.File(OutFile, 'w')
    
    grp = fOut.create_group("Extrapolated_N2.dir")
    
    ## Compute for only (2,2) modes
    if only22:
      print("Computing for (2,2) mode only")

      modes = [(2,2)]
      for mode in modes:
        l = mode[0]
        m = mode[1]
        print("Computing for ", mode)
        time, total = ComputedCSModifiedStrain(p, mode, ell, lev)
        ## Compute for the given mode
        dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
        (len(time),3), dtype='f')

        dataset[:,0] = time
        dataset[:,1] = np.real(total)
        dataset[:,2] = np.imag(total)

      ## For all other modes
      l_arr = range(2, 9)
      for l in l_arr:
        print("Computing for l = ", l)

        for m in range(-l, l+1):
          mode = (l, m)
          if mode not in modes:
            print("Setting zero for ", mode)
            dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
            (len(time),3), dtype='f')

            dataset[:,0] = np.zeros_like(time)
            dataset[:,1] = np.zeros_like(time)
            dataset[:,2] = np.zeros_like(time)
            
    ## All modes except m = 0
    elif dropm0:
        print("Computing for all modes except m = 0")
    
        l_arr = range(2, 9)

        for l in l_arr:
            print("Computing for l = ", l)

            for m in range(-l, l+1):
            
                mode = (l, m)
        
                if m != 0:
                    ## Compute the mode if not m = 0
                    time, total = ComputedCSModifiedStrain(p, mode, ell, lev)

                    dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
                    (len(time),3), dtype='f')

                    dataset[:,0] = time
                    dataset[:,1] = np.real(total)
                    dataset[:,2] = np.imag(total)
            
                else:
                    ## Set to zero if m = 0 
                    print("Setting zero for ", mode)
                    dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
                    (len(time),3), dtype='f')

                    dataset[:,0] = np.zeros_like(time)
                    dataset[:,1] = np.zeros_like(time)
                    dataset[:,2] = np.zeros_like(time)
        
    ## Compute for all of the modes
    else:
        print("Computing for all of the modes")
        l_arr = range(2, 9)

        for l in l_arr:
            print("Computing for l = ", l)

            for m in range(-l, l+1):
                mode = (l, m)
                print(mode)

                ## Compute for the given mode
                time, total = ComputedCSModifiedStrain(p, mode, ell, lev)

                dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
                (len(time),3), dtype='f')

                dataset[:,0] = time
                dataset[:,1] = np.real(total)
                dataset[:,2] = np.imag(total)

    fOut.close()
    print("Wrote waveforms to file", OutFile)

def GetModesFromString(modes):
    """ Method to get output modes from a given string, used to specify
        the modes to the LVC file generation """
    if modes == 'all':
        modes = [[l,m] for l in range(2,9) for m in range(-l,l+1)]
    elif modes == '22only':
        modes = [[2, 2]]
    return modes

def GenerateStrainFiles(ell, lev, only22, dropm0, convert_lvc):
    """ Generates the sxs format waveform and LVC format waveform for a given
        beyond-GR simulation with coupling parameter ell. 
        
        The file BeyondGRAnalysis/Waveforms/Lev*/BackgroundStrain.h5 contains the strain
        of the GR background, h_GR. 
        
        Meanwhile the file BeyondGRAnalysis/Waveforms/Lev*/DeltaStrain.h5 contains the 
        leading-order modification to the beyond-GR strain, delta_h. 
        
        We first generate the sxs format total beyond-GR waveform for these
        using 
        
        h = h_GR + ell^4 * delta_h. 
        
        Then we convert the resulting waveform into the LVC format, from 
        which we can then generate frames files
        
        only22 is a boolean specifying if we should only output the 22 mode, 
        which can be useful for testing purposes
        
        """
    ## Since ell has a decimal point and *.h5 readers cannot
    ## handle this decimal point, we will replace the decimal point 
    ## with the character p, ie
    ## 0.1 will become 0p1
    ell_string = str(ell).replace('.', 'p')

    ## Call to generate total waveform in sxs format
    OutputdCSModifiedStrain("Waveforms/", ell, only22, dropm0, lev)

    if convert_lvc:
      ## Vars needed for the sxs to lvc conversion script. 
      ## values like resolution and sxs_metadata and sxs_resolutions
      ## don't really matter, but we need to specify them
      sxs_data = "Waveforms/Lev" + str(lev) ## input directory
      resolution = 0
      sxs_metadata = "catalog_tools/Metadata/sxs_catalog.json"
      sxs_resolutions = "catalog_tools/Metadata/sxs_catalog_resolutions.json"
      ## all modes are included in the file structure -- even though we may have
      ## set all except (2,2) or all of the m = 0 modes to zero
      modes = GetModesFromString("all") ## modes to output
      out_path = "Waveforms/dCS_" + ell_string + "_Lev" + str(lev) ## output directory
      in_file = "dCS_" + ell_string + "_Lev" + str(lev) + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
      out_file = "Waveforms/dCS_" + ell_string + "_Lev" + str(lev) + "/dCS_ell_" + ell_string + ".h5"

      ## convert the simulation into sxs format
      convert_simulation(sxs_data, resolution, sxs_metadata, sxs_resolutions, modes, out_path, \
                    in_name = in_file, out_name = out_file)
      print("Output LVC format waveform to", out_file)
    
def main():
  p = argparse.ArgumentParser(description="Generate dCS waveform for a given coupling parameter value")
  p.add_argument("--ell", required=True, type=float,\
    help="Value of dCS coupling constant")
  p.add_argument("--lev", required=True, type=int,\
    help="Resolution of NR run to use")
  p.add_argument('--only22', help='Only output the 22 mode', \
    dest='only22', action='store_true')
  p.add_argument('--dropm0', help='Include all modes up to l = 8 except m = 0 modes', \
    dest='dropm0', action='store_true')
  p.add_argument('--convert_lvc', help='Output waveform in LVC format as well', \
    dest='convert_lvc', action='store_true')
  p.set_defaults(only22=False)
  p.set_defaults(dropm0=False)
  p.set_defaults(convert_lvc=False)
  args = p.parse_args()

  GenerateStrainFiles(args.ell, args.lev, args.only22, args.dropm0, args.convert_lvc)

if __name__ == "__main__":
  main()
