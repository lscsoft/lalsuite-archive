#!/bin/bash

################################################################################
# edit these appropriately

month_gps_time='847555570'
month_duration='2419200'
cat='CAT_3'

coire_path='/home/cdcapano/local/s5_2yr_lowcbc_20080829/bin/lalapps_coire'
plotifar_path='/home/cdcapano/local/cbc_s5_1yr_20070129/bin/plotifar'
data_type='playground_only'

log_path='/usr1/cdcapano/log'
condor_priority='20'

# don't touch anything below here
################################################################################

#gps_end_time is needed for plotifar (don't need to edit)
gps_end_time=$(( ${month_gps_time} + ${month_duration} ))

#generate dag
/bin/echo -n "Generating ifar_result_${data_type}.dag and .sub files..."
if [ 1 ]; then
  #apply data filter coire
  for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    #to double time files
    for combo in H1L1 H2L1; do
      glob_files="corse_all_data_files/${mass}/${combo}-CORSE_ALL_DATA_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      outfile="ifar_result_files/${data_type}/${mass}/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      summaryfile="ifar_result_files/summ_files/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.txt"
      echo "JOB $outfile ifar_result_${data_type}.filter_coire.sub"
      echo "RETRY $outfile 1"
      echo "VARS $outfile macroglob=\"$glob_files\" macrooutfile=\"$outfile\" macrosummaryfile=\"$summaryfile\""
      echo "CATEGORY $outfile filter_coire"
      echo
    done
    #to triple time files
    for combo in H1L1 H2L1 H1H2L1; do
      glob_files="corse_all_data_files/${mass}/H1H2L1_${combo}-CORSE_ALL_DATA_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      outfile="ifar_result_files/${data_type}/${mass}/separated_triple_time/H1H2L1_${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      summaryfile="ifar_result_files/summ_files/H1H2L1-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.txt"
      echo "JOB $outfile ifar_result_${data_type}.filter_coire.sub"
      echo "RETRY $outfile 1"
      echo "VARS $outfile macroglob=\"$glob_files\" macrooutfile=\"$outfile\" macrosummaryfile=\"$summaryfile\""
      echo "CATEGORY $outfile coire"
      echo
    done
    #coire together triple time files
    glob_files="ifar_result_files/${data_type}/${mass}/separated_triple_time/*-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outfile="ifar_result_files/${data_type}/${mass}/H1H2L1-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    echo "JOB $outfile ifar_result_${data_type}.coire.sub"
    echo "RETRY $outfile 1"
    echo "VARS $outfile macroglob=\"$glob_files\" macrooutfile=\"$outfile\""
    echo "CATEGORY $outfile coire"
    for combo in H1L1 H2L1 H1H2L1; do
      parent_file="ifar_result_files/${data_type}/${mass}/separated_triple_time/H1H2L1_${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      echo "PARENT $parent_file CHILD $outfile"
    done
    echo
    #run plotifar on the individual mass bins
    for combo in H1L1 H2L1 H1H2L1; do
      job_name="${combo}-plotifar_${mass}_${data_type}"
      glob_file="ifar_result_files/${data_type}/${mass}/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      outpath="ifar_result_files/${data_type}/${mass}"
      if [[ "$combo" == "H1H2L1" ]]; then
        triggers="--h1-triggers --h2-triggers --l1-triggers"
        num_cat='3' #because 3 types of triggers in triple time
      elif [[ "$combo" == "H1L1" ]]; then
        triggers="--h1-triggers --l1-triggers"
        num_cat='1'
      elif [[ "$combo" == "H2L1" ]]; then
        triggers="--h2-triggers --l1-triggers"
        num_cat='1'
      fi
      user_tag="${mass}-${data_type}"
      t_correct_num="ifar_result_files/summ_files/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.txt"
      t_correct_denom="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${combo}-${month_gps_time}-${month_duration}.txt"
      echo "JOB $job_name ifar_result_${data_type}.plotifar.sub"
      echo "RETRY $job_name 1"
      echo "VARS $job_name macroglob=\"$glob_file\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotriggers=\"$triggers\" macronumcat=\"$num_cat\"  macrousertag=\"$user_tag\" macrotcorrnum=\"$t_correct_num\" macrotcorrdenom=\"$t_correct_denom\""
      echo "CATEGORY $job_name plotifar"
      parent_file="ifar_result_files/${data_type}/${mass}/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      echo "PARENT $parent_file CHILD $job_name"
      echo
    done
  done
  #coire together mass files
  for combo in H1L1 H2L1 H1H2L1; do
    glob_files="ifar_result_files/${data_type}/mchirp_*/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outfile="ifar_result_files/${data_type}/${combo}-CORSE_ALL_MASSES_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    echo "JOB $outfile ifar_result_${data_type}.coire.sub"
    echo "RETRY $outfile 1"
    echo "VARS $outfile macroglob=\"$glob_files\" macrooutfile=\"$outfile\""
    echo "CATEGORY $outfile coire"
    for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
      parent_file="ifar_result_files/${data_type}/${mass}/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      echo "PARENT $parent_file CHILD $outfile"
    done
    echo
    #run plotifar on CORSE_ALL_MASSES files
    job_name="${combo}-plotifar_${data_type}"
    glob_file="ifar_result_files/${data_type}/${combo}-CORSE_ALL_MASSES_${data_type}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outpath="ifar_result_files/${data_type}"
    if [[ "$combo" == "H1H2L1" ]]; then
      triggers="--h1-triggers --h2-triggers --l1-triggers"
      num_cat='9' #because 3 types of triggers in triple time * 3 mass bins
    elif [[ "$combo" == "H1L1" ]]; then
      triggers="--h1-triggers --l1-triggers"
      num_cat='3' #because 1 type of trigger in double time * 3 mass bins
    elif [[ "$combo" == "H2L1" ]]; then
      triggers="--h2-triggers --l1-triggers"
      num_cat='3' #because 1 type of trigger in double time * 3 mass bins
    fi
    t_correct_num="ifar_result_files/summ_files/${combo}-CORSE_${data_type}_${cat}-${month_gps_time}-${month_duration}.txt"
    t_correct_denom="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${combo}-${month_gps_time}-${month_duration}.txt"
    echo "JOB $job_name ifar_result_${data_type}.plotifar.sub"
    echo "RETRY $job_name 1"
    echo "VARS $job_name macroglob=\"$glob_file\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotriggers=\"$triggers\" macronumcat=\"$num_cat\"  macrousertag=\"$data_type\" macrotcorrnum=\"$t_correct_num\" macrotcorrdenom=\"$t_correct_denom\""
    echo "CATEGORY $job_name plotifar"
    echo "PARENT $outfile CHILD $job_name"
    echo
  done
  echo "MAXJOBS filter_coire 20"
  echo "MAXJOBS coire 20"
  echo "MAXJOBS plotifar 20"
fi > ifar_result_${data_type}.dag

if [ 1 ]; then
  echo "universe = standard"
  echo "executable = ${coire_path}"
  echo "arguments = --glob \$(macroglob) --output \$(macrooutfile) --data-type ${data_type} --summary \$(macrosummaryfile)"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/coire-\$(cluster)-\$(process).err"
  echo "output = logs/coire-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > ifar_result_${data_type}.filter_coire.sub

if [ 1 ]; then
  echo "universe = standard"
  echo "executable = ${coire_path}"
  echo "arguments = --glob \$(macroglob) --output \$(macrooutfile) --data-type ${data_type}"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/coire-\$(cluster)-\$(process).err"
  echo "output = logs/coire-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > ifar_result_${data_type}.coire.sub

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${plotifar_path}"
  echo "arguments = --glob \$(macroglob) --output-path \$(macrooutpath) --ifo-times \$(macroifotimes) \$(macrotriggers) --gps-start-time ${month_gps_time} --gps-end-time ${gps_end_time} --num-categories \$(macronumcat) --combine-types --ifar-dist --enable-output --user-tag \$(macrousertag) --t-correct-num \$(macrotcorrnum) --t-correct-denom \$(macrotcorrdenom)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/plotifar-\$(cluster)-\$(process).err"
  echo "output = logs/plotifar-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > ifar_result_${data_type}.plotifar.sub

#make directory structure
if [ ! -d ifar_result_files ] ; then
 mkdir ifar_result_files
fi
if [ ! -d ifar_result_files/summ_files ] ; then
  mkdir ifar_result_files/summ_files
fi
if [ ! -d ifar_result_files/${data_type} ] ; then
  mkdir ifar_result_files/${data_type}
fi
for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
  if [ ! -d ifar_result_files/${data_type}/${mass} ] ; then
    mkdir ifar_result_files/${data_type}/${mass}
  fi
  if [ ! -d ifar_result_files/${data_type}/${mass}/separated_triple_time ] ; then
    mkdir ifar_result_files/${data_type}/${mass}/separated_triple_time
  fi
done
echo " done."
echo "*************************************************************"
echo "  Now run: condor_submit_dag ifar_result_${data_type}.dag"
echo "*************************************************************"

