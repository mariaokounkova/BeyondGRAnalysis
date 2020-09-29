## Take in the value of ell as a command-line argument
## SNR is the second command-line argument
## The third command-line argument is the NR resolution of the run
## replace . in the string with p in order
## to be able to read from the .h5 files
ell=$1
ell=`echo ${ell/./'p'}`
dist=$2
lev=$3

DIR="dCS_"${ell}"_Dist_"${dist}"_Lev"${lev}
h5path=Waveforms_22/${DIR}/dCS_ell_${ell}.h5

NAME="dCS_"${ell}

PSDs=PSDs/design

# Move files into directory
mv H-H1HWINJ_${NAME}.gwf Waveforms_22/${DIR}/H-H1HWINJ_${NAME}.gwf
mv L-L1HWINJ_${NAME}.gwf Waveforms_22/${DIR}/L-L1HWINJ_${NAME}.gwf

# make the cache files
echo -e "-\t-\t-\t-\tfile://localhost/home/maria.okounkova/BeyondGRAnalysis/Waveforms_22/${DIR}/H-H1HWINJ_${NAME}.gwf" > Waveforms_22/${DIR}/H-${NAME}.lcf
echo -e "-\t-\t-\t-\tfile://localhost/home/maria.okounkova/BeyondGRAnalysis/Waveforms_22/${DIR}/L-L1HWINJ_${NAME}.gwf" > Waveforms_22/${DIR}/L-${NAME}.lcf

# now make the Bayeswave ini file
cp dCS.ini Waveforms_22/${DIR}/${NAME}.ini
sed -i "s/ELL/${ell}/g" "Waveforms_22/${DIR}/${NAME}.ini"
sed -i "s/NAME/${NAME}/g" "Waveforms_22/${DIR}/${NAME}.ini"
sed -i "s/LEV/${lev}/g" "Waveforms_22/${DIR}/${NAME}.ini"

# make the Bayeswave submission script
cp run_bw.sh Waveforms_22/${DIR}/run_bw_${NAME}.sh
sed -i "s/NAME/${NAME}/g" "Waveforms_22/${DIR}/run_bw_${NAME}.sh"
