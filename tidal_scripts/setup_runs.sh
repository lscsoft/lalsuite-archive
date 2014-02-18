# Load your tidal_dev environment
. /home/magathos/lalsuites/LoadEnv_tidal_dev.sh

# Define seed, run and log directory (on LHO, LLO, CIT, log should be /usr1/... instead of /home/... ; on Atlaas it is /local/user/...)
seed=5000
RUN_DIR=/home/magathos/Tidal/runs/TF2_15_lso/$seed
LOG_DIR=/usr1/magathos/Tidal/runs/TF2_15_lso/$seed

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

for inj in {MS1,H4,SQM3,PP}
do 
  lalapps_inspinj --f-lower 20.0 --gps-start-time 932172000 --gps-end-time 933172000 --seed ${seed} --waveform TaylorF2threePointFivePN --min-distance 1.00e+05 --max-distance 2.50e+05 --d-distr volume --l-distr random --i-distr uniform --m-distr componentMass --min-mass1 1.0 --max-mass1 2.0 --min-mass2 1.0 --max-mass2 2.0 --min-mtotal 2.0 --max-mtotal 4.0 --time-step 1.0e+03 --amp-order 0 --disable-spin --eos ${inj} --output inj_${inj}_${seed}.xml
  mkdir inj_${inj}
  cd inj_${inj}
  for rec in {MS1,H4,SQM3,PP}
  do 
    mkdir rec_${rec}
    cd rec_${rec}
    lalinference_pipe -I ${RUN_DIR}/inj_${inj}_${seed}.xml -p ${LOG_DIR}/inj_${inj}/rec_${rec} -r ${RUN_DIR}/inj_${inj}/rec_${rec} ${RUN_DIR}/parser_${rec}.ini
    cd ..
  done
  cd ..
done
