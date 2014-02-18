# Load your tidal_dev environment
. /home/vdbroeck/env_tidal_dev.sh

# Define seed, run and log directory (on LHO, LLO, CIT, log should be /usr1/... instead of /home/... ; on Atlaas it is /local/user/...)
seed=5000
RUN_DIR=/home/vdbroeck/Tidal_contact_spins_gaussian_tidal_dev/$seed
LOG_DIR=/usr1/vdbroeck/Tidal_contact_spins_gaussian_tidal_dev/$seed

if [ ! -d ${RUN_DIR} ]
then
  mkdir -p ${RUN_DIR}
fi

if [ ! -d ${LOG_DIR} ]
then
  mkdir -p ${LOG_DIR}
fi

cp ./parser_*.ini ${RUN_DIR}/

cd ${RUN_DIR}

for inj in MS1 H4 SQM3 PP
do 
  lalapps_inspinj --f-lower 40.0 --gps-start-time 932172000 --gps-end-time 933172000 --seed ${seed} --waveform TaylorF2threePointFivePN --min-distance 1.00e+05 --max-distance 2.50e+05 --d-distr volume --l-distr random --i-distr uniform --m-distr gaussian --min-mass1 1.0 --max-mass1 1.7 --min-mass2 1.0 --max-mass2 1.7 --min-mtotal 2.0 --max-mtotal 3.4 --mean-mass1 1.35 --stdev-mass1 0.05 --mean-mass2 1.35 --stdev-mass2 0.05 --time-step 1.0e+03 --amp-order 0 --enable-spin --aligned --spin-gaussian --mean-spin1 0.0 --stdev-spin1 0.05 --mean-spin2 0.0 --stdev-spin2 0.05 --min-spin1 0.0 --max-spin1 0.7 --min-spin2 0.0 --max-spin2 0.7 --eos ${inj} --output inj_${inj}_${seed}.xml
  mkdir inj_${inj}
  cd inj_${inj}
  for rec in MS1 H4 SQM3 PP
  do 
    mkdir rec_${rec}
    cd rec_${rec}
    lalinference_pipe -I ${RUN_DIR}/inj_${inj}_${seed}.xml -p ${LOG_DIR}/inj_${inj}/rec_${rec} -r ${RUN_DIR}/inj_${inj}/rec_${rec} ${RUN_DIR}/parser_${rec}.ini
    cd ..
  done
  cd ..
done
