import numpy as np
import argparse
import math

def get_wavelet_params(filename, model, chirpflag=False, O1version=False, **keyword_parameters):
    """
    Read in chain file and get all wavelet params
    
    arguments
    ---------
    filename (str): the chain file
    
    model (str): signal or glitch
    
    optional, chirpflag: True if using chirplets
    optional, O1version: True if using O1 era chains
    optional, restrict (int): line number if you only want one draw from the chain
    
    outputs
    -------
    dictionary of the wavelet params
    """
    NW = 5 # number of intrinsic parameters (changes for chirplets)
    NE = 6 # number of extrinsic parameters
    start = 1
    
    labels = ['t','f','Q','A','phi']
    extlabels = ['alpha','sindelta','psi','ecc', 'phi0','scale']
    if chirpflag:
        NW = 6
        labels.append('beta')
    
    data = {}
    for l in labels:
        data[l] = []
    
    if model == 'signal': # get intrinsic parameters
        for l in extlabels:
            data[l] = []

    data['D'] = []

    infile = open(filename)
    lines = infile.readlines()

    if ('restrict' in keyword_parameters):
        restrict = int(keyword_parameters['restrict'])
        rn = [restrict]
    else:
        rn = np.arange(0,len(lines))


    for j in rn:
        line = lines[j]
        spl = line.split()
        waveletnumber = int(spl[0]) # how many wavelets
        data['D'].append(waveletnumber)
        if model == 'signal':
            start = NE+1 # extra parameters
            if O1version:
                start += 1
            for l in range(0,NE):
                data[extlabels[l]].append(float(spl[l+1]))

        for i in range(0,waveletnumber):
            for l in range(0,NW):
                data[labels[l]].append(float(spl[start+i*NW+l]))

    return data

def wt(wave_params,psdfile):
    """
    Makes a waveform from a set of wavelets
    
    arguments
    ---------
    wave_params (dict): the wavelet parameters
    
    psdfile (str): data file of the PSD
    
    outputs
    -------
    array of the waveform
    """
    psd = np.genfromtxt(psdfile)
    
    # Find what Nsamp should be:
    l = len(psd)
    l = 2*l
    Nsamp = 2**(l-1).bit_length()
    
    hs = np.zeros(Nsamp)
    
    # Find Tobs (1/df)
    Tobs = 1./(psd[1,0]-psd[0,0])
    
    fmin = int(psd[0,0]*Tobs)
    
    
    wavenumber = wave_params['D'][0]
    
    for j in range(0,wavenumber):
    
        t0 = wave_params['t'][j]
        f0 = wave_params['f'][j]
        Q = wave_params['Q'][j]
        A = wave_params['A'][j]
        phi0 = wave_params['phi'][j]
        
        i = int(f0*Tobs)
        fac = 1.0/math.sqrt(psd[i-fmin,1])
        
        tau = Q/(2*np.pi*f0)
        
        tmax = t0 + 4.0*tau
        tmin = t0 - 4.0*tau
        
        imin = int((tmin/Tobs)*(Nsamp))
        imax = int((tmax/Tobs)*(Nsamp))
        if imin < 0: imin = 0
        if imax > Nsamp: imax = Nsamp
        
        for i in range(imin,imax):
            t = float(i)/Nsamp*Tobs
            sf = A*np.exp(-((t-t0)**2)/(tau**2))
            sf *= fac
            hs[i] += sf*np.cos(2*np.pi*f0*(t-t0)+phi0)

    return hs

################################
# Deal with chirplets later
################################

#def wt_chirplets(hs,wave_params,psd):
#    
#    wavenumber = wave_params['D']
#    
#    for j in range(0,wavenumber):
#        
#        t0 = wave_params['t'][j]
#        f0 = wave_params['f'][j]
#        Q = wave_params['Q'][j]
#        A = wave_params['A'][j]
#        phi0 = wave_params['phi'][j]
#        
#        i = int(f0*Tobs)
#        fac = 1.0/math.sqrt(psd[i])
#        
#        tau = Q/(2*np.pi*f0)
#        
#        tmax = t0 + 4.0*tau
#        tmin = t0 - 4.0*tau
#        
#        imin = int((tmin/Tobs)*(NMAX))
#        imax = int((tmax/Tobs)*l(NMAX))
#        if imin < 0: imin = 0
#        if imax > NMAX: imax = NMAX
#        
#        for i in range(imin,imax):
#            t = i/NMAX*Tobs
#            sf = A*np.exp(-((t-t0)**2)/(tau**2))
#            sf *= fac
#            hs[i] += sf*np.cos(2*np.pi*f0*(t-t0)+phi0)



