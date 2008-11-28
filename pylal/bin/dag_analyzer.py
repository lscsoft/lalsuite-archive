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
#################################
class DAGParser:

  # ------------------------------------------------
  def __init__(self, filename):

    # Storing filename and the file's content
    self.filename = filename
    self.content = None
    
    # job dict and status
    self.job_dict = {}
    self.status = -1

    # lists to store the node numbers
    self.vector_node_time = []
    self.vector_node_number = []

    # variables to store the run status and start/end times
    self.starttime = []
    self.endtime = []
    self.wallclocktime = 0
    self.actual_idle_jobs = 0
    self.queued_jobs = 0
    self.running_jobs = 0
    self.number_done = 0
    self.number_idle = 0
    self.number_ready = 0
    self.number_unready = 0
    self.number_failed = 0
    
    self.verbose = False

    self.def_status = {-1:"Unknown",0:"Unstarted",1:"Running", 2:"ERROR",\
                       3:"FINISHED"}

  # ------------------------------------------------
  def add_job(self, id):
    self.job_dict[id]={'id':id, 'start':-1, 'end':None, \
                  'startVec':[],'endVec':[],'sub': None, 'time':0}

  # ------------------------------------------------
  def convert_date(self, date_string, time_string):

    month = int(date_string.split('/')[0])
    day = int(date_string.split('/')[1])
    string = "08/%02d/%02d %s" % (month, day, time_string)
    
    time_class = time.strptime(string, "%y/%m/%d %H:%M:%S")
    time_sec = time.mktime( time_class )
    
    return time_sec


  # ------------------------------------------------
  def get_status(self):
    
    try:
      f = open(self.filename)
      self.content = f.readlines()
      f.close()
      self.status = 1  # file exists

      if self.verbose:
        print "Parsing file ", self.filename

      # parsing file
      self.parse()
      
    except IOError:
      self.status = 0 # file does not exist: DAG has not started

    if self.status == 1:
      if 'EXITING' in self.content[-1]:
        if 'STATUS 0' in self.content[-1]:
          self.status = 3 # DAG stopped: Finish
        if 'STATUS 1' in self.content[-1]:
          self.status = 2 # DAG stopped: Error

    
  # ------------------------------------------------
  def parse(self):

    # get the wall-time start/end times of this outfile
    #words = self.content[0].split()
    #self.starttime = self.convert_date(words[0], words[1])
    #words = self.content[-1].split()
    #self.endtime = self.convert_date(words[0], words[1])

    # reset some counters
    actual_node_number = 0
    flag_counter = -1

    # and parse the file, line for line...
    for line in self.content:  
      words = line.split()

      # create the real wallclock-running time
      if 'STARTING UP' in line:
        self.starttime.append(self.convert_date(words[0], words[1]))
      if 'EXITING' in line:
        self.endtime.append(self.convert_date(words[0], words[1]))          

      ## start of a job
      if 'ULOG_EXECUTE' in line:

        time_sec = self.convert_date(words[0], words[1])
        id = words[7]

        # add a new job to the dictionary if not already contained
        if not self.job_dict.has_key(id):
          self.add_job(id)

        # store the maximum start time (because jobs might get evicted)
        if self.job_dict[id]['start'] is not None:
          self.job_dict[id]['start'] = max(  self.job_dict[id]['start'], time_sec)
        self.job_dict[id]['startVec'].append( time_sec )

      ## end of a job
      if 'ULOG_JOB_TERMINATED' in line:

        time_sec = self.convert_date(words[0], words[1])
        id = words[7]

        # store it
        if not self.job_dict.has_key( id):
          self.add_job(id)
        self.job_dict[id]['end'] = time_sec
        self.job_dict[id]['endVec'].append( time_sec )

      ## end of a job
      if 'ULOG_JOB_EVICTED' in line:

        time_sec = self.convert_date(words[0], words[1])
        id = words[7]

        # store it if new
        if not self.job_dict.has_key( id):
          print "WARNING: Job with ID %s is being evicted, but has never started before!"\
                % id
          self.add_job(id)
        self.job_dict[id]['endVec'].append(time_sec)

      ## type of the job
      if 'submitting' in line:

        for word in words:
          if 'dag_node_name' in word:
            index = words.index(word)

        # get ID and the submit-type
        id =  words[index+2].split("'")[1]
        sub = words[len(words)-1]

        # store it
        if not self.job_dict.has_key( id):
          self.add_job(id)
        self.job_dict[id]['sub']=sub


      ## get the number of RUNNING jobs  
      flag_counter -= 1

      if 'Number of idle job procs' in line:
        self.actual_idle_jobs = int(words[7])

      if 'Done' in line:
        flag_counter = 2

      if flag_counter==0:

        # get the time and the number of running jobs in this line
        time_sec = self.convert_date(words[0], words[1])
        self.queued_jobs = int(words[4])    
        self.running_jobs = self.queued_jobs - self.actual_idle_jobs

        self.number_done = int(words[2])
        self.number_idle = self.actual_idle_jobs
        self.number_ready =  int(words[6])
        self.number_unready = int(words[7])
        self.number_failed = int(words[8])        

        # and store it in the list
        self.vector_node_time.append(time_sec)
        self.vector_node_number.append(self.running_jobs)


    # check the wallclock time
    if len(self.starttime) != len(self.endtime):
      if self.status != 1:
        print "Some kind of error with the time"

      words = self.content[-1].split()
      self.endtime.append(self.convert_date(words[0], words[1]))

    for start, end in zip(self.starttime, self.endtime):
      self.wallclocktime += end-start

      

  def print_status(self):
    """
    Puts together the status of this DAG. 
    """
    
    text = self.def_status[self.status]

    # see if DAG is running, add run informations
    if self.status == 1:

      text+=" | Done: %d  Executed: %d  idle: %d  Unready: %d  Failed: %d " %\
             (self.number_done, self.running_jobs, self.number_idle+self.number_ready,\
              self.number_unready, self.number_failed)
      
    return text

  
#################################  
def check_subdags(dag_file, directory, verbose = False):
  
  file  = open(dag_file)
  content = file.readlines()
  file.close()

  # one DAG we know about ...
  filenames = [dag_file]

  for line in content:
    if "JOB" in line and "DIR" in line:
      strs = line.split()

      # set the directory
      if not "DONE" in line:
        fname = directory+'/'+strs[-1] + "/"
      else:
        print "Seems someone tweaked this file! This sub-DAG"\
              "will not be considered."
        continue
      
      for string in strs:
        if "condor.sub" in string:
          fname += string[:-11]

      # Append to list
      filenames.append(fname)
      if verbose:
        print 'Found sub-DAG ' + fname

  return filenames

  
#################################  


## make the glob
directory = './'
if len(sys.argv)>1:
  user_input = sys.argv[1]
  directory = os.path.dirname(user_input)
  if user_input.endswith('dag'):
    filelist = glob.glob(user_input)
  else:
    filelist = glob.glob(user_input+'*.dag')    
else:
  filelist = glob.glob('*dag')


## check if there is only one DAG (or Uber-DAG)
dag_file = None
if len(filelist)>1:
  print "WARNING: Two dags found. Use the one with 'uberdag' in it...\n"
  for file in filelist:
    if 'uberdag' in file:
      dag_file = file
else:
  dag_file = filelist[0]

## check for sub-dags
dag_list = check_subdags(dag_file, directory)


## check the status of the DAG's itself
verbose_all = False
create_plots = False


## MASTER loop over all found DAG's and sub-DAG's
list_dag = []
all_dag_time = 0.0
c = 0
for dag_name in dag_list:

  c += 1
  out_name = dag_name + '.dagman.out'

  # parse the DAG in the DAG Parser class
  dag = DAGParser(out_name)
  dag.get_status()

  ## final analysis of all the jobs from ALL out files

  # check if there are any 'problematic' jobs in the DAG
  number_correct = 0
  number_incorrect = 0

  ## retrieve the complete running time
  for id, job in dag.job_dict.iteritems():

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
      #print "The following job has problems: (WHAT PROBLEMS??)"
      #print job


  ## loop over the different job-types:
  total_time = 0.0
  total_jobs = 0
  for job_type in ['inspiral_H1','inspiral_H2','inspiral_L1',\
                   'inspiral_veto_H1','inspiral_veto_H2','inspiral_veto_L1',\
                   'trigbank','tmpltbank',\
                   'coire','thinca_','thinca2']:

    done_job = False
    list_duration = []
    for id, job in dag.job_dict.iteritems():
      if job['sub'] is None:
        continue
      if job['end'] and job_type in job['sub']:
        done_job = True
        list_duration.append(job['time'])

    if done_job:
      if verbose_all:
        print "%18s:  %4d jobs run with a mean duration of %6.0f seconds." % \
              (job_type, len(list_duration), numpy.mean(list_duration) )
      total_time += sum(list_duration)
      total_jobs += len(list_duration)

      if create_plots:
        data =  numpy.asarray(list_duration)
        pylab.clf()
        pylab.hist(data/3600.0)
        pylab.grid(True)
        pylab.xlabel('Time of the job to run [h]',size='x-large')
        pylab.ylabel('number',size='x-large')
        pylab.savefig( 'dag_analyzer-hist_'+job_type+'.png')


  ## master print out
  if verbose_all:
    print "\nThe total number of jobs found is %d" % total_jobs
    print "There are %d jobs for which no information is available" % number_incorrect

  ## get the runing times for all DAG's 
  runtime = dag.wallclocktime
  all_dag_time += runtime
  print "\n----------------------------------------------------------"
  print "DAG %s: %s" % (os.path.basename(dag_name), dag.print_status())
  print "Total analysis time: %.1f CPU hours (%.1f CPU days) on ~%.1f nodes| "\
        "Wallclock: %.1f hours (%.1f days)" % \
        ( total_time/3600.0, total_time/86400.0, numpy.mean(dag.vector_node_number), \
          runtime/3600.0, runtime/86400.0)
  if verbose_all:
    print "-- The total number of jobs found is %d" % total_jobs
    print "-- There are %d jobs for which no information is available" % number_incorrect

  if create_plots:
    # create a plot of the used nodes during the run(s)
    px = numpy.asarray(vector_node_time)-dag_starttime[0]
    py = numpy.asarray(vector_node_number)
    pylab.clf()
    pylab.plot( px/86400.0,  py, 'r.')
    pylab.grid(True)
    pylab.xlabel('Days since start of DAG',size='x-large')
    pylab.xlabel('Number of computer nodes used',size='x-large')
    pylab.savefig('dag_analyzer-usedNodes.png')

print "\n=========================================================="
print "The TOTAL running time of all DAG's is %.1f hours (%.1f days)" % \
      (all_dag_time/3600.0, all_dag_time/86400.0)
