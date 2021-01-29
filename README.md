**Using and Generating Beyond-GR waveforms**

To generate the BackgroundStrain.h5 and DeltaStrain.h5 files, run commands like

python3 Process_dCS_Waveforms.py --waveform_dir /home/maria.okounkova/BeyondGRAnalysis/Waveforms/Lev2

To then generate dCS strain waveforms for a certain resolution and coupling constant, do

python3 Generate_dCS_Strain.py --ell 0 --lev 2

To then make the frames files, use the notebook

GenerateFramesFiles.ipynb

To then run Bayeswave, go into the Frames directory, and do 

./run_bw.sh

