#!/bin/bash

################################################################################
# get needed options from ini file

data_type=`cat write_ifar_scripts.ini | grep 'data_type' | awk '{print $3}'`

month_gps_time=`cat write_ifar_scripts.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts.ini | grep 'cat' | awk '{print $3}'`

coire_path=`cat write_ifar_scripts.ini | grep 'coire_path' | awk '{print $3}'`
combinedFAR_path=`cat write_ifar_scripts.ini | grep 'combinedFAR_path' | awk '{print $3}'`

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
  # run plotifar on all the mass bins globbed together 
  for combo in H1L1 H2L1 H1H2L1; do
    job_name="${combo}-ifar-combined_ALL_MASSES_${data_type}_slide"
    glob_files="corse_all_data_files/${data_type}_slide/${combo}*CORSE_SLIDE_`echo ${data_type} | tr '[a-z]' '[A-Z]'`_*_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outpath="combined_ifar_files/${data_type}_slide/"
    time_correct_file="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
    user_tag="ALL_MASSES-${data_type}_slide"
    echo "JOB $job_name ifar_result.ifar_combined.sub"
    echo "RETRY $job_name 0"
    echo "VARS $job_name macroglob=\"$glob_files\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotcorrfile=\"$time_correct_file\" macrousertag=\"$user_tag\""
    echo "CATEGORY $job_name ifarcombined"
    echo
  done
  echo "MAXJOBS ifarcombined 20"
fi > combined_ifar_${data_type}_slide.dag

if [ 1 ]; then
  # run plotifar on all the mass bins globbed together 
  for combo in H1L1 H2L1 H1H2L1; do
    job_name="${combo}-ifar-combined_ALL_MASSES_${data_type}"
    glob_files="corse_all_data_files/${data_type}/${combo}*CORSE_`echo ${data_type} | tr '[a-z]' '[A-Z]'`_*_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outpath="combined_ifar_files/${data_type}/"
    time_correct_file="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
    user_tag="ALL_MASSES-${data_type}"
    echo "JOB $job_name ifar_result.ifar_combined.sub"
    echo "RETRY $job_name 0"
    echo "VARS $job_name macroglob=\"$glob_files\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotcorrfile=\"$time_correct_file\" macrousertag=\"$user_tag\""
    echo "CATEGORY $job_name ifarcombined"
    echo
  done
  echo "MAXJOBS ifarcombined 20"
fi > combined_ifar_${data_type}.dag


if [ 1 ]; then
  for injstring in BBHINJ BNSINJ NSBHINJ SPININJ; do
    # run plotifar on all the mass bins globbed together 
    for combo in H1L1 H2L1 H1H2L1; do
      job_name="${combo}-ifar-combined_ALL_MASSES_${injstring}"
      glob_files="corse_all_data_files/${injstring}/${combo}*CORSE_${injstring}_*_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      outpath="combined_ifar_files/${injstring}/"
      time_correct_file="second_coire_files/summ_files_${injstring}/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
      user_tag="ALL_MASSES-${injstring}"
      echo "JOB $job_name ifar_result.ifar_combined.sub"
      echo "RETRY $job_name 0"
      echo "VARS $job_name macroglob=\"$glob_files\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotcorrfile=\"$time_correct_file\" macrousertag=\"$user_tag\""
      echo "CATEGORY $job_name plotifar"
      echo
    done
  done
  echo "MAXJOBS plotifar 20"
fi > combined_ifar_injection.dag

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${combinedFAR_path}"
  echo "arguments = --glob \$(macroglob) --output-path \$(macrooutpath) --ifo-times \$(macroifotimes) --gps-start-time ${month_gps_time} --gps-end-time ${gps_end_time} --time-correct-file \$(macrotcorrfile) --user-tag \$(macrousertag)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/ifarcombined-\$(cluster)-\$(process).err"
  echo "output = logs/ifarcombined-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > ifar_result.ifar_combined.sub

#make directory structure
if [ ! -d combined_ifar_files ] ; then
 mkdir combined_ifar_files
fi
if [ ! -d combined_ifar_files/summ_files ] ; then
  mkdir combined_ifar_files/summ_files
fi
for string in ${data_type} ${data_type}_slide BNSINJ BBHINJ NSBHINJ SPININJ; do
  if [ ! -d combined_ifar_files/${string} ] ; then
    mkdir combined_ifar_files/${string}
  fi
done
echo " done."
echo "*******************************************************************"
echo "  Now run: condor_submit_dag combined_ifar_${data_type}.dag"
echo "      and: condor_submit_dag combined_ifar_${data_type}_slide.dag"
echo "      and: condor_submit_dag combined_ifar_injection.dag"
echo "*******************************************************************"

