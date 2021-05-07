#!/usr/bin/env python
"""
Tutorial to demonstrate running parameter estimation on a reduced parameter
space for an injected signal.

This example estimates the masses using a uniform prior in both component masses
and distance using a uniform in comoving volume prior on luminosity distance
between luminosity distances of 100Mpc and 5Gpc, the cosmology is Planck15.
"""

import numpy as np
import bilby

# Set the duration and sampling frequency of the data segment that we're
# going to inject the signal into
duration = 2.
sampling_frequency = 2048.

# Specify the output directory and the name of the simulation.
outdir = '/home/maria.okounkova/BeyondGRAnalysis/BilbyPE/Surrogate/outdir'
label = 'fast_tutorial'
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.

# Going through all of the parameter meanings
# a_1, a_2: spin magnitude 0 <= a <= 1
# tilt_1, tilt_2: 0 for straight up, pi for straight down (co-latitude to spin vector in orbital frame)
# theta_jn: not the inclination of orbit to line of sight, but inlicnation of total angular momentum
# to the line of sight. (Farr and Farr have a paper) cos(theta_jn) = total angular momentum (j) . vector pointing to observer (n)
# if 0, then pointing at you. pi, looking at the bottom. 
# psi: polarization (angle on the sky of the pericenter wrt to interferometer -- degenerate for circular orbits)
# phase: phase of coalescence (at some fiducial time)
# phi_jl: longitudes, azimuthal angles. Angle between total angular momentum projected into the x-y plane. 0 for non-precession --
# should maybe be called phi_j (azimuthal angle of total angular momentum vector in orbital frame)
# phi_12: Angle from spin 1 to spin 2 also in x-y plane in orbital frame 
# Normal to the plane and orbital angular momentum are not the same (bc of PN corrections). May not be the instantaneous plane where
# the velocities of the two objects live. 

## Injection in detector frame quantities, and recovers in detector frame quantities 

## Flip the mass ratio
injection_parameters = dict(
    mass_ratio = 0.8, chirp_mass = 29.378924880880138, a_1=0.33, a_2=0.44, 
    tilt_1=0.0, tilt_2=3.14159265359,
    phi_12=0.0, phi_jl=0.0, luminosity_distance=2000., theta_jn=3.14159265359, psi=0.0,
    phase=0.0, geocent_time=1126259462.0, ra=1.952318922, dec=-1.26967171703)



# Fixed arguments passed into the source model
#waveform_arguments = dict(waveform_approximant='IMRPhenomPv2',
#                          reference_frequency=50., minimum_frequency=20.)
waveform_arguments = dict(waveform_approximant='NRSur7dq4',
                          reference_frequency=25., minimum_frequency=25.)
print("Set up waveform arguments")

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    #time_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments)
print("Set up waveform generator")

# Set up interferometers.  In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 0.5)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)
print("Injected signal")

# Set up a PriorDict, which inherits from dict.
# By default we will sample all terms in the signal models.  However, this will
# take a long time for the calculation, so for this example we will set almost
# all of the priors to be equal to their injected values.  This implies the
# prior is a delta function at the true, injected value.  In reality, the
# sampler implementation is smart enough to not sample any parameter that has
# a delta-function prior.
# The above list does *not* include mass_1, mass_2, theta_jn and luminosity
# distance, which means those are the parameters that will be included in the
# sampler.  If we do nothing, then the default priors get used.
priors = bilby.gw.prior.BBHPriorDict()

priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=injection_parameters['geocent_time'] - 1,
    maximum=injection_parameters['geocent_time'] + 1,
    name='geocent_time', latex_label='$t_c$', unit='$s$')
#for key in ['luminosity_distance', 'theta_jn', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'psi', 'ra',
#            'dec', 'geocent_time', 'phase']:
for key in ['luminosity_distance', 'theta_jn', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'psi', 'ra',
            'dec', 'geocent_time', 'phase']:
    priors[key] = injection_parameters[key]
# If we're searching for the mass ratio, restrict the prior to the allowed values for NRSur7dq4
priors['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.2, maximum=1, latex_label='$q$')
# If we're searching for the spins, restrict ourselves to 0.8
priors['a_1'] = bilby.core.prior.Uniform(name='a_1', minimum=0.0, maximum=0.7, latex_label='$a_1$')
priors['a_2'] = bilby.core.prior.Uniform(name='a_2', minimum=0.0, maximum=0.7, latex_label='$a_2$')

# Initialise the likelihood by passing in the interferometer data (ifos) and
# the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator)

# Run sampler.  In this case we're going to use the `dynesty` sampler
result = bilby.run_sampler(
    likelihood=likelihood, priors=priors, sampler='dynesty', npoints=1000,
    injection_parameters=injection_parameters, outdir=outdir, label=label)

# Make a corner plot.
result.plot_corner()

