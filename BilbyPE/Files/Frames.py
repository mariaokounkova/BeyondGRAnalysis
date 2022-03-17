#!/usr/bin/env python

import bilby
from gwpy.timeseries import TimeSeries
from bilby.gw import utils as gwutils
import numpy as np

######################
#### Setup ###########
######################

root_dir = '/home/maria.okounkova/BeyondGRAnalysis/BilbyPE/FRAMES_DIR/'
logger = bilby.core.utils.logger
outdir = root_dir + 'outdir'
label = 'fast_tutorial' 
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# Interferometer injections
interferometer_names = ['H1', 'L1','V1']
ifo_list = bilby.gw.detector.InterferometerList([])

######################
#### Data Reading ####
######################

geocenter_time = 1126259462.0
duration = 4  # length of data segment containing the signal
time_of_event = 1126259462.0
post_trigger_duration = 2  # time between trigger time and end of segment
end_time = time_of_event + post_trigger_duration
start_time = end_time - duration

for det in interferometer_names:

    logger.info("Downloading analysis data for ifo {}".format(det))
    file_name = root_dir + "Frames/" + det + ".gwf"
    ifo = bilby.gw.detector.get_empty_interferometer(det)

    ## Read in the strain data from frame
    data = gwutils.read_frame_file(file_name = file_name, \
                            start_time = start_time, end_time = end_time, \
                            buffer_time = 4.0, \
                            channel = det + ":LDAS_STRAIN")

    ifo.set_strain_data_from_gwpy_timeseries(data)
    
    ## Read in our own psd
    
    psd_file = root_dir + 'aLIGOZeroDetHighPower-PSD_25Hz.txt'
    psd_frequencies_value, psd_value = np.loadtxt(psd_file, comments="#",usecols=([0,1]),unpack=True)

    ## Comment out if you want to use the default design psd
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=psd_frequencies_value, psd_array=psd_value)
    
    ifo_list.append(ifo)


logger.info("Saving data plots to {}".format(outdir))
bilby.core.utils.check_directory_exists_and_if_not_mkdir(outdir)
ifo_list.plot_data(outdir=outdir, label=label)

######################
#### Priors ##########
######################

##### CHOOSE PRIOR FILE
##### DEFAULT PRIOR FILES: GW150914.prior, binary_black_holes.prior,
##### Needs to specify path if you want to use any other prior file.
prior = bilby.gw.prior.BBHPriorDict(filename = root_dir +'Frame.prior')
logger.info("Set up priors")

######################
#### Analysis ########
######################

sampling_frequency = 2048.  # same at which the data is stored 

# Set up waveform arguments
# NRSur7dq4: Surrogate model
# SEOBNRv4PHM: EOB waveform from GWTC-3 analysis
# IMRPhenomXPHM: Phenom waveform from GWTC-3 analysis
waveform_arguments = dict(waveform_approximant='IMRPhenomXPHM',
                          minimum_frequency = 25.0,
                          reference_frequency = 30.0,
                          maximum_frequency = 2048.0)
logger.info("Set up waveform arguments")

# Waveform generator 
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    #time_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments)
logger.info("Set up waveform generator")

# Likelihood
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifo_list, waveform_generator=waveform_generator)
logger.info("Set up likelihood")

# Run sampler 
## TODO: Get injection parameters in as arg injection_parameters=injection_parameters
result = bilby.run_sampler(
    likelihood=likelihood, priors=prior, sampler='dynesty', npoints=1080,
    outdir=outdir, label=label, npool=18)
logger.info("Ran sampler")

result.plot_corner()
