#!/usr/bin/env python
"""
"""
import bilby
from gwpy.timeseries import TimeSeries
from bilby.gw import utils as gwutils
import numpy as np

######################
#### Setup ###########
######################

logger = bilby.core.utils.logger
outdir = '/home/maria.okounkova/BeyondGRAnalysis/BilbyPE/Surrogate_Frames/outdir'
label = 'fast_tutorial' 
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# Interferometer injections
interferometer_names = ['H1', 'L1']
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
print(start_time, end_time)

roll_off = 0.4  # smoothness in a tukey window, default is 0.4s
# This determines the time window used to fetch open data
psd_duration = 1 * duration
psd_start_time = start_time - psd_duration
psd_end_time = start_time
print(psd_start_time, psd_end_time)
filter_freq = None  

for det in interferometer_names:

    logger.info("Downloading analysis data for ifo {}".format(det))
    file_name = "/home/maria.okounkova/BeyondGRAnalysis/BilbyPE/" + \
                "Surrogate_Frames/Frames_EqualMassNonSpinning/" \
                + det + ".gwf"
    ifo = bilby.gw.detector.get_empty_interferometer(det)

    ## Read in the strain data from frame
    data = gwutils.read_frame_file(file_name = file_name, \
                            start_time = start_time, end_time = end_time, \
                            buffer_time = 4.0, \
                            channel = det + ":LDAS_STRAIN")

    ifo.set_strain_data_from_gwpy_timeseries(data)

    logger.info("Downloading psd data for ifo {}".format(det))
    ## Read in the fake psd data (all zeros in our case)
    psd_data = gwutils.read_frame_file(file_name = \
                    "/home/maria.okounkova/BeyondGRAnalysis/BilbyPE/" + \
                    "Surrogate_Frames/Frames_PSD/" \
                                    + det + ".gwf", \
                          start_time = psd_start_time, \
                          end_time = psd_end_time, \
                          buffer_time = 4.0, \
                          channel = det + ":LDAS_STRAIN")


    psd_alpha = 2 * roll_off / duration  # shape parameter of tukey window
    psd = psd_data.psd(
        fftlength=duration,
        overlap=0,
        window=("tukey", psd_alpha),
        method="median")
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=psd.frequencies.value, psd_array=psd.value)
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
prior = bilby.gw.prior.BBHPriorDict(filename='Frame.prior')

######################
#### Analysis ########
######################

sampling_frequency = 2048.  # same at which the data is stored ## Masha changed from 4096

# Set up waveform arguments
waveform_arguments = dict(waveform_approximant='NRSur7dq4',
                          reference_frequency=25., minimum_frequency=25.)
#waveform_arguments = {
#    'waveform_approximant': 'IMRPhenomPv2',
#    'reference_frequency': 50  # most sensitive frequency
#}

# Waveform generator 
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    #time_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments)

#waveform_generator = bilby.gw.WaveformGenerator(
#    parameter_conversion=conversion,
#    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
#    waveform_arguments=waveform_arguments)

# Likelihood
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifo_list, waveform_generator=waveform_generator)

#likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
#    interferometers=ifo_list, waveform_generator=waveform_generator,
#    priors=prior, time_marginalization=False, distance_marginalization=False,
#    phase_marginalization=False)

# Run sampler 
## TODO: Get injection parameters in as arg injection_parameters=injection_parameters
result = bilby.run_sampler(
    likelihood=likelihood, priors=prior, sampler='dynesty', npoints=1000,
    outdir=outdir, label=label)

# Implemented Samplers:
# LIST OF AVAILABLE SAMPLERS: Run -> bilby.sampler.implemented_samplers
#npoints = 512  # number of live points for the nested sampler
#n_steps = 100  # min number of steps before proposing a new live point,
# defaults `ndim * 10`
#sampler = 'dynesty'
# Different samplers can have different additional kwargs,
# visit https://lscsoft.docs.ligo.org/bilby/samplers.html for details.

# result = bilby.run_sampler(
#     likelihood, prior, outdir=outdir, label=label,
#     sampler=sampler, nlive=npoints, use_ratio=False,
#     walks=n_steps, n_check_point=10000, check_point_plot=True,
#     conversion_function=bilby.gw.conversion.generate_all_bbh_parameters)

result.plot_corner()
