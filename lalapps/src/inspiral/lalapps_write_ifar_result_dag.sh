#!/bin/bash

################################################################################
# get needed options from ini file

data_type=`cat write_ifar_scripts.ini | grep 'data_type' | awk '{print $3}'`

month_gps_time=`cat write_ifar_scripts.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts.ini | grep 'cat' | awk '{print $3}'`

coire_path=`cat write_ifar_scripts.ini | grep 'coire_path' | awk '{print $3}'`
plotifar_path=`cat write_ifar_scripts.ini | grep 'plotifar_path' | awk '{print $3}'`

log_path=`cat write_ifar_scripts.ini | grep 'log_path' | awk '{print $3}'`
condor_priority=`cat write_ifar_scripts.ini | grep 'condor_priority' | awk '{print $3}'`

#Print options out to screen for verification
echo "Options used are:"
echo "  data_type = ${data_type}"
echo "  month_gps_time = ${month_gps_time}"
echo "  month_duration = ${month_duration}"
echo "  cat = ${cat}"
echo "  coire_path = ${coire_path}"
echo "  plotifar_path = ${plotifar_path}"
echo "  log_path = ${log_path}"
echo "  condor_priority = ${condor_priority}"
echo
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
      t_correct_denom="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
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
    t_correct_denom="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
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

