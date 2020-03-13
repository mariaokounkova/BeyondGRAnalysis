#!/usr/bin/python

import argparse
import scipy
import h5py
import numpy as np
import scipy.integrate
import scipy.interpolate
import json
import sxs
from scipy.optimize import fmin

def CutTimes(time, data, TLow, TUp): 
	""" Cut time and data to be between 
	    TLow and TUp  """
	TLowIndex = np.where(time >= TLow)[0][0]
	TUpIndex = np.where(time <= TUp)[0][-1]
	time = time[TLowIndex:TUpIndex]
	data = data[TLowIndex:TUpIndex]
	return time, data

def GetPeakTimeMode(time, data): 
	""" Grab the peak time of some data """
	t_peak = time[np.argmax(np.absolute(data))]
	return t_peak

def SubtractPeakTimeMode(time, data): 
	""" Subtract the peak time of some data """
	t_peak = GetPeakTimeMode(time, data)
	return time - t_peak

def InterpolateTimes(time, data, time_dest):
    """ Interpolates time, data onto new time axis
        time_dest """
    ## build the interpolant only in the region where we need it
    ## time, data = CutTimes(time, data, min(time_dest), max(time_dest))
    interpolant = scipy.interpolate.CubicSpline(time, data)
    return interpolant(time_dest)

def RunningIntegral(time, data):
    """ Given an array, compute an integral of the array between
        the start time and each point in time"""
    ## Go over each time and compute the integral to that time
    integral = np.zeros(len(time))
    for n in range(1, len(time)):
        t = time[0:n]
        d = data[0:n]
        integ = scipy.integrate.simps(d, x=t)
        integral[n] = integ
    return time, integral

def DeltaPsi4Factor(hPsi4, B = 0.1):
    """ h^{(2)}_{ab} = (ell/GM)^4 / 8 Delta g_{ab} 
        -- need to divide by a factor of 8
    """
    hPsi4 = hPsi4 / 8.0
    ## Also need to divide by B^2 (since it's at order l^4)
    hPsi4 = hPsi4 / B**2
    return hPsi4

def ReadExtrapolatedMode(p, piece, mode, order=2, ell=None):
	""" Given a file of extrapolated modes, read in the (mode)
	    at a given order """
	ell = name = str(ell).replace('.', 'p')
	piece_dict = {"hRWZ" : "/rhOverM_Asymptotic_GeometricUnits.h5", \
				  "DeltaStrain" : "/DeltaStrain.h5", \
				  "BackgroundStrain" : "/BackgroundStrain.h5", \
				  "Psi4" : "/rMPsi4_Asymptotic_GeometricUnits.h5", \
				  "DeltaPsi4" : "/rMDeltaPsi4_Asymptotic_GeometricUnits.h5", \
				  "dCSModified" : "/rhOverM_Asymptotic_GeometricUnits_dCS_ell_" + ell + ".h5"}
	try:
		file = p + piece_dict[piece]
	except: file = p + piece
	l = mode[0]
	m = mode[1]
	f = h5py.File(file, 'r')
	data = f['Extrapolated_N'+str(order)+'.dir']['Y_l' + str(l) + '_m'  + str(m) + '.dat']
	time, re, im = data[:,0], data[:,1], data[:,2]
	result = re + 1j*im
	## Put in the factor if it's DeltaPsi4
	if piece == "DeltaPsi4":
		result = DeltaPsi4Factor(result)
	return time, result

def ChristodolouMass(data_dir):
	""" Grab the sum of the Christodolou masses from 
	    metadata.json file """
	with open(data_dir + 'metadata.json') as json_file:
		data = json.load(json_file)
		m1 = data['reference_mass1']
		m2 = data['reference_mass2']
		return m1 + m2

def ComputeStrain(time, psi4, data_dir):
    """ Compute the delta strain from a psi4 -- 
        could be background psi4, or delta psi4, or 
        background psi4 + l**4 * delta psi4. 
        Integrate both the real and imag parts 
        Have to do each part separately because scipy integrate doesn't 
        understand integrating a complex number. Luckily integration is a 
        linear operation """
    ## Grab the sum of the Christodolou masses

    mass = ChristodolouMass(data_dir) 
    ## Divide the result twice by the Christodolou mass, since we 
    ## have rDeltaMPsi4, and we need rhOverM
    psi4_real = np.real(psi4)
    psi4_imag = np.imag(psi4)

    news_real_time, news_real = RunningIntegral(time, psi4_real)
    news_imag_time, news_imag = RunningIntegral(time, psi4_imag)

    strain_real_time, strain_real = RunningIntegral(news_real_time, news_real)
    strain_imag_time, strain_imag = RunningIntegral(news_imag_time, news_imag)

    strain = (strain_real + 1j*strain_imag)/ mass**2
    return strain_real_time, strain

def GetStrainFromPsi4(p, mode, order=2):
	""" Generate the background strain from background psi4 file. 
	    Useful if instead of the background RWZ strain we want to 
	    use the integrated psi4 in order to be consistent with how 
	    we're computing delta strain"""
	time, psi4 = ReadExtrapolatedMode(p, "Psi4", mode, order)
	strain_time, strain = ComputeStrain(time, psi4, p)
	return strain_time, strain

def ComputedCSDeltaStrain(p, mode):
	""" Given a path p to the extrapolated hPsi4, 
	    compute the modification to the gravitational wave strain. 
	    Return time time, strain, and delta_strain """

	print("Computing delta strain for " + str(mode) + "\n")
	## Grab delta psi4
	delta_time, delta_psi4 = ReadExtrapolatedMode(p, "DeltaPsi4", mode)

	## Grab the background strain
	time, strain = ReadExtrapolatedMode(p, "hRWZ", mode)

	## Compute delta strain
	delta_time, delta_strain = ComputeStrain(delta_time, delta_psi4, p)
	
	## First check that we do in fact start delta_strain at zero
	if(delta_strain[0] != 0):
		print("Waveform perturbation is initially non-zero!")
		return

	## Pad the delta strain array to be the length of the simmulation
	delta_start_time = delta_time[0]

	## Cut the time array
	time_before = time[np.where(time < delta_start_time)]
	time_after  = time[np.where(time >= delta_start_time)]

	## Make an array of zeroes for the time before
	delta_strain_pad = np.zeros(len(time_before))

	## Now add in the cut array to the delta array
	delta_time = np.concatenate((time_before, delta_time))
	delta_strain = np.concatenate((delta_strain_pad, delta_strain))

	## Now working with the padded array
	## Interpolate onto the perturbed time axis
	#strain = InterpolateTimes(time, strain, delta_time)
	delta_strain = InterpolateTimes(delta_time, delta_strain, time)

	return time, strain, delta_strain
	
def OutputdCSDeltaStrain(p, only22):
	""" Given a path p to the extrapolated hPsi4, 
	    compute the modification to the gravitational wave strain. 
	    Return time time, strain, and delta_strain """

	## Dump the result to a file
	OutFile1 = p + 'BackgroundStrain.h5'
	OutFile2 = p + 'DeltaStrain.h5'

	fOut1 = h5py.File(OutFile1, 'w')
	fOut2 = h5py.File(OutFile2, 'w')
    
	grp1 = fOut1.create_group("Extrapolated_N2.dir")
	grp2 = fOut2.create_group("Extrapolated_N2.dir")
    
	l_arr = range(2, 9) if not only22 else [2]

	for l in l_arr:
		print("Computing for l =", l)
		for m in range(-l, l+1) if not only22 else [2, -2]:

			mode = (l, m)

			## Compute strain and delta strain for this mode
			time, strain, delta_strain = ComputedCSDeltaStrain(p, mode)

			dataset1 = grp1.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
				(len(time),3), dtype='f')
			dataset2 = grp2.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
				(len(time),3), dtype='f')
			dataset1[:,0] = time
			dataset1[:,1] = np.real(strain)
			dataset1[:,2] = np.imag(strain)

			dataset2[:,0] = time
			dataset2[:,1] = np.real(delta_strain)
			dataset2[:,2] = np.imag(delta_strain)

	return time, strain, delta_strain

def ComputedCSModifiedStrain(p, mode, l):
	""" Given a value of the dCS coupling constant l, a path 
	    p to the extrapolated hPsi4 compute the modified gravitational wave strain """

	## Read in the background
	time, strain = ReadExtrapolatedMode(p, "BackgroundStrain", mode)
	delta_time, delta_strain = ReadExtrapolatedMode(p, "DeltaStrain", mode)

	## Now add the strain and delta strain together
	## with the correct value of l
	total = strain + l**4 * delta_strain

	print("Computing strain for l = ", l)
	return time, total


def OutputdCSModifiedStrain(p, ell, only22):
    """ Generate an h5 file with the modified strain for a 
        given value of ell """

    ## For naming the file, replace . with p because otherwise
    ## the .h5 file can't be read by catalog scripts
    name = str(ell).replace('.', 'p')
    
    OutFile = p + 'rhOverM_Asymptotic_GeometricUnits_dCS_ell_' + name + '.h5'
    fOut = h5py.File(OutFile, 'w')
    
    grp = fOut.create_group("Extrapolated_N2.dir")
    
    l_arr = range(2, 9) if not only22 else [2]

    for l in l_arr:
    	print("Computing for l =", l)
    	for m in range(-l, l+1) if not only22 else [2, -2]:

    		mode = (l, m)

    		## Compute for the given mode
    		time, total = ComputedCSModifiedStrain(p, mode, ell)

    		dataset = grp.create_dataset("Y_l"+str(l)+"_m"+str(m)+".dat", \
				(len(time),3), dtype='f')

    		dataset[:,0] = time
    		dataset[:,1] = np.real(total)
    		dataset[:,2] = np.imag(total)

    fOut.close()

def Overlap(data1, data2): 
	""" Compute overlap between two waveforms in 
	    a certain time window. This is done 
	    accroding to the conventions in the SXS catalog
	    paper, including taking the real part of the 
	    complex mismatch """
	def product(d1, d2):
		return np.dot(d1, np.conj(d2))

	denom = np.sqrt(product(data1, data1) * product(data2, data2))
	numerator = product(data1, data2)
	mismatch = 1.0 - np.real(numerator/denom)
	return mismatch

def OverlapInRegion(time1, data1, time2, data2, t_min, t_max, length=10000):
	""" Compute the overlap between two waveforms in a certain time region, 
	    by interpolating the data onto a time axis linearly spaced between
	    t_min and t_max with length length """
	time_window = np.linspace(t_min, t_max, length)
	data1 = InterpolateTimes(time1, data1, time_window)
	data2 = InterpolateTimes(time2, data2, time_window)
	mismatch = Overlap(data1, data2)
	return mismatch

def ComputeMinOverlap(t_min, t_max, time1, data1, time2, data2):
	""" Compute the overlap minimizing over time shift """
	def f(shift):
		## compute the overlap for this data
		overlap = OverlapInRegion(time1 + shift, data1, time2, data2, t_min, t_max)
		return overlap

	minimum = fmin(f, x0=0)
	return f(minimum)

def ReadLVCStrainMode(p, ell, mode, order=2): 
    """ Read in the the phase and amplitude of the strain from 
        the LVC format file """
    l = mode[0]
    m = mode[1]
    h_file = p + "dCS_ell_" + \
                 str(ell).replace('.', 'p') + ".h5"
    print("Reading in the LVC format strain from:", h_file)
    f = h5py.File(h_file, 'r')
    amp_dat = f["amp_l"+str(l)+"_m"+str(m)]
    phase_dat = f["phase_l"+str(l)+"_m"+str(m)]
    time = amp_dat['X']
    amp = amp_dat['Y']
    phase = phase_dat['X']
    return time, amp, phase
    
def MakeJSonFile(data_dir):
	a = sxs.metadata.Metadata()
	b = a.from_txt_file(data_dir + 'metadata.txt')

def main():
	p = argparse.ArgumentParser(description="Generate dCS waveform from given dCS simulation")
	#p.add_argument("--ell", required=True, type=float,\
	#	help="Value of dCS coupling constant")
	p.add_argument("--waveform_dir", required=True, \
		help="Directory containing extrapolated, snipped simulation waveforms.")
	p.add_argument('--only22', help='Only output the 22 mode', \
		dest='only22', action='store_true')
	p.set_defaults(only22=False)
	args = p.parse_args()

	data_dir = str(args.waveform_dir) + '/'

	print("\n Working in directory " + data_dir + "\n")

	## First make the json file
	print("Making metadata.json")
	MakeJSonFile(data_dir + '/')

	## produce DeltaStrain.h5 and BackgroundStrain.h5 files
	print("\n Computing and outputting delta strain \n")
	OutputdCSDeltaStrain(data_dir, args.only22)


if __name__ == "__main__":
  main()
