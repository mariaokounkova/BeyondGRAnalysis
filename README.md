Using and Generating Beyond-GR waveforms 

All of the instructions below will correspond to the directory

/home/maria.okounkova/BeyondGR

on CIT

All dynamical Chern-Simons gravity waveforms are in 

/home/maria.okounkova/BeyondGR/Waveforms/dCS_*

where the * corresponds to the value of the dCS coupling constant, 
with '.' replaced by 'p' (so dCS_0p1 corresponds to dCS_0.1). Note that
dCS_0p0 corresponds to the GR waveform.

Data is presently available for dCS coupling constant values

0.0 (GR), 0.1, 0226 (Max allowed by perturbative scheme)

and SNRs (for aligo design sensitivity)

20, 50, 80, 100, 150

But I include instructions below for how to make your own combinations
if you want

----------------------------------------------

(1) If you want SXS format strain waveforms (up to modes l = 8), 
these are available in each waveform directory as

Waveforms/dCS_*/rhOverM_Asymptotic_GeometricUnits_dCS_ell_*.h5

where again * corresponds to the value of the dCS coupling constant.
These are in the standard SXS format, and thus can be processed with
any scripts in https://github.com/sxs-collaboration/catalog_tools, 
for example.

----------------------------------------------

(2) If you want LVC format strain waveforms, these are available in 
each waveform directory as

waveforms/dCS_*/dCS_ell_*.h5

This intermediate format is needed to go between the SXS format
and the frames files. 

----------------------------------------------

(3) If you want to generate your own dCS waveform with a 
desired coupling constant, simply do

python3 Generate_dCS_Strain.py --ell [value of dCS coupling constant]

(also see python3 Generate_dCS_Strain.py -h for help)

Note that this takes a while because a spline interpolant must be
build for each mode.

This will generate a directory

Waveforms/dCS_* with both the SXS format and LVC format waveforms. 

----------------------------------------------

(4) If you want frames files, these are in each waveform directory as

Waveforms/dCS_*/L-L1HWINJ_dCS_*_SNR_**.gwf
Waveforms/dCS_*/H-H1HWINJ_dCS_*_SNR_**.gwf

Where * corresponds to the value of the dCS coupling constant, 
and ** corresponds to the SNR. 

----------------------------------------------

(5) If you want to generate your own frames files with a desired
value of the dCS coupling constant and SNR, run 

./create_frame_file_from_NR_data.sh [ell] [SNR]

where the command-line arguments are the desired dCS 
coupling constant and the desired SNR.

Note that the LVC format waveforms for the desired ell must
already exist in Waveforms/dCS_* (see (3) for instructions 
on how to generate this)

Note that this will also create a Bayeswave .ini file
and Bayeswave condor submission script

----------------------------------------------

(6) If you want to run Bayeswave, the ini file is in each
waveform directory as

Waveforms/dCS_*/dCS_*_SNR_**.ini

where again * is the desired value of the dCS coupling constant
and ** is the SNR

The corresponding script to submit a condor job to run Bayeswave is

Waveforms/dCS_*/run_bw_dCS_*_SNR_**.sh

----------------------------------------------

(7) If you want to make a Bayeswave run for a desired dCS coupling constant
and SNR, then follow steps (3), (5), and (6)


