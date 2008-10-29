#!/usr/bin/python

import os
import sys
import glob
import time

import matplotlib
matplotlib.use('Agg')
import pylab 
import numpy



#################################
def convert_date(date_string, time_string):

  month = int(date_string.split('/')[0])
  day = int(date_string.split('/')[1])
  string = "08/%02d/%02d %s" % (month, day, time_string)
  
  time_class = time.strptime(string, "%y/%m/%d %H:%M:%S")
  time_sec = time.mktime( time_class )

  return time_sec

#################################  
def add_job(job_dict, id):
  job_dict[id]={'id':id, 'start':-1, 'end':None, \
                'startVec':[],'endVec':[],'sub': None, 'time':0}
  
      
#################################
def parse_outfile(job_dict, all_lines):

  # get the wall-time start/end times of this outfile
  words = all_lines[0].split()
  dag_starttime.append(convert_date(words[0], words[1]))
  words = all_lines[-1].split()
  dag_endtime.append(convert_date(words[0], words[1]))

  # reset some counters
  actual_node_number = 0
  flag_counter = -1

  # and parse the file, line for line...
  for line in all_lines:  
    words = line.split()
    
    ## start of a job
    if 'ULOG_EXECUTE' in line:

      time_sec = convert_date(words[0], words[1])
      id = words[7]

      # add a new job to the dictionary if not already contained
      if not job_dict.has_key(id):
        add_job( job_dict, id )

      # store the maximum start time (because jobs might get evicted)
      if job_dict[id]['start'] is not None:
        job_dict[id]['start'] = max(  job_dict[id]['start'], time_sec)
      job_dict[id]['startVec'].append( time_sec )

    ## end of a job
    if 'ULOG_JOB_TERMINATED' in line:

      time_sec = convert_date(words[0], words[1])
      id = words[7]

      # store it
      if not job_dict.has_key( id):
        add_job( job_ict, id )
      job_dict[id]['end']=time_sec
      job_dict[id]['endVec'].append( time_sec )

    ## end of a job
    if 'ULOG_JOB_EVICTED' in line:

      time_sec = convert_date(words[0], words[1])
      id = words[7]

      # store it if new
      if not job_dict.has_key( id):
        print ">>> This line should never be called !!!"
        add_job( job_dict, id )
      job_dict[id]['endVec'].append( time_sec )

    ## type of the job
    if 'submitting' in line:

      for word in words:
        if 'dag_node_name' in word:
          index=  words.index(word)

      # get ID and the submit-type
      id =  words[index+2].split("'")[1]
      sub = words[len(words)-1]

      # store it
      if not job_dict.has_key( id):
        add_job( job_dict, id )
      job_dict[id]['sub']=sub


    ## get the number of RUNNING jobs  
    flag_counter -= 1
    
    if 'Number of idle job procs' in line:
      actual_idle_jobs = int(words[7])

    if 'Done' in line:
      flag_counter = 2
    
    if flag_counter==0:

      # get the time and the number of running jobs in this line
      time_sec = convert_date(words[0], words[1])
      queued_jobs = int(words[4])    
      running_jobs = queued_jobs - actual_idle_jobs
      
      # and store it in the list
      vector_node_time.append(time_sec)
      vector_node_number.append(running_jobs)



  
  
#################################  


## make the glob
if len(sys.argv)>1:
  directory = sys.argv[1]
else:
  directory = './'

filelist = glob.glob(directory+'*.dagman.out')

# init the vectors and the data structures
vector_node_time = []
vector_node_number = []
job_dict = dict()
dag_starttime = []
dag_endtime = []

## loop over all dagman.out files
for filename in filelist:

  # try to open the file and read the lines
  try:

    # open file and read all lines
    file = open(filename)
    all_lines =  file.readlines()
    file.close()

    # parse this file
    print "Parsing file ", filename
    parse_outfile(job_dict, all_lines)
    
  except IOError:
    raise "File not found: ", filename



## final analysis of all the jobs from ALL out files

number_correct = 0
number_incorrect = 0

## retrieve the complete running time
for id, job in job_dict.iteritems():

  if job['end']:
    number_correct += 1
    time_end =  job['end']
    time_start =  job['start']        
    duration = time_end - time_start

    if len( job['startVec'] )>1:
      # job got evicted; adding the single run times
      # TODO: Is this the right way? Yes, its the total running time,
      # even if job starts from scratch
      duration = 0
      for start, end in zip( job['startVec'],  job['endVec']):
        duration +=end-start

    job['time'] = duration
    
  else:
    number_incorrect += 1
    #print "The following job has problems:"
    #print job


print "\n"

## loop over the different job-types:
total_time = 0.0
total_jobs = 0
for job_type in ['inspiral_H1','inspiral_H2','inspiral_L1',\
                 'inspiral_veto_H1','inspiral_veto_H2','inspiral_veto_L1',\
                 'trigbank','tmpltbank',\
                 'coire','thinca_','thinca2']:

  done_job = False
  list_duration = []
  for id, job in job_dict.iteritems():
    if job['end'] and job_type in job['sub']:
      done_job = True
      list_duration.append(job['time'])

  if done_job:
    print "%18s:  %4d jobs run with a mean duration of %6.0f seconds." % \
          (job_type, len(list_duration), numpy.mean(list_duration) )
    total_time += sum(list_duration)
    total_jobs += len(list_duration)
    

    data =  numpy.asarray(list_duration)
    pylab.clf()
    pylab.hist(data/3600.0)
    pylab.grid(True)
    pylab.xlabel('Time of the job to run [h]',size='x-large')
    pylab.ylabel('number',size='x-large')
    pylab.savefig( 'dag_analyzer-hist_'+job_type+'.png')

print "\nThe total number of jobs found is %d" % total_jobs
print "There are %d jobs for which no information is available" % number_incorrect
print "The total analysis time: %.1f CPU hours (%.1f CPU days)  " % \
      (total_time/3600.0, total_time/86400.0)
print "\n"

## analyze the total dag times
total_dag_time = 0.0
for i in range(len(dag_starttime)):
  starttime = dag_starttime[i]
  endtime = dag_endtime[i]
  runtime = endtime-starttime
  total_dag_time += runtime
  
  print "The %dth DAG took %.1f hours (%.1f days) to run " % \
        (i+1, runtime/3600.0, runtime/86400.0)
  
print "** The TOTAL running time of all DAG's is %.1f hours (%.1f days)" % \
      (total_dag_time/3600.0, total_dag_time/86400.0)


# create a plot of the used nodes during the run(s)
px = numpy.asarray(vector_node_time)-dag_starttime[0]
py = numpy.asarray(vector_node_number)
pylab.clf()
pylab.plot( px/86400.0,  py, 'r.')
pylab.grid(True)
pylab.xlabel('Days since start of DAG',size='x-large')
pylab.xlabel('Number of computer nodes used',size='x-large')
pylab.savefig('dag_analyzer-usedNodes.png')
