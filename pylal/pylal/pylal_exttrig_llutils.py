#!/usr/bin/python

import os
import sys
import pickle
import time
import subprocess
import ConfigParser
import optparse
import pickle
import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')

from glue import segments
from glue import segmentsUtils
from glue.ligolw import lsctables
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import utils
from pylal.plotsegments import PlotSegmentsPlot
from pylal.grbsummary import multi_ifo_compute_offsource_segment as micos
from pylal import antenna

# the config parser to be used in some of the functions
cp = None

template_trigger_hipe = "./lalapps_trigger_hipe"\
  " --h1-segments H1-science.txt" \
  " --l1-segments L1-science.txt" \
  " --v1-segments V1-science.txt" \
  " --list %s" \
  " --grb %s --onsource-left 5" \
  " --onsource-right 1 --config-file %s" \
  " --injection-config injectionsWI.ini --log-path %s" \
  " --number-buffer-left 8 --number-buffer-right 8" \
  " --num-trials 340 --padding-time 72 --verbose" \
  " --skip-datafind --skip-dataquality"

ifo_list = ['H1','L1','V1']

offset_gps_to_linux = 315982785

total_summary_prefix = """
<body style="color: rgb(0, 0, 0); background-color: rgb(221, 255, 255);" alink="#000099" link="#000099" vlink="#990099">

<h1>Summary of Gamma Ray Burst low-latency results during S6</h1>

<span style="font-weight: bold;"><br><br>
The following table contain a list of Gamma Ray Bursts occured during S6, with information about time, position on the sky, as well as duration and redshift (if available). This table has been automatically created by pylal_exttrig_llmonitor (in pylal_exttrig_llutils.py) to show a summary of the low-latency inspiral analysis of the GRBs during S6. A page describing this search can be found in the <a href="https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/S6Plan/090706044855TriggeredSearchLow_Latency_Exttrig_Search#preview">wiki</a>.<br><br>

The number in the IFO columns indicate the antenna factor for this GRB and this detector, if there is data available at this time (1 meaning optimal location, 0 meaning worst location). In addition, a <b>bold</b> number indicates that the current segment is long enough to contain the required number of off-source segments around the GRB, not required to be symmetrical.<br><br>
Date of last creation: %s<br><br>

</span><span style="font-weight: bold;">
<br><br>
</div>
<table border="1" cellpadding="2" cellspacing="2">
  <tbody>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Nr</td>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">GRB</td>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Status</td>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">GPS<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Date<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">redshift<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">duration [s]<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">RA<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">DEC<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">H1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">L1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">V1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Result<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Link<br>
"""

# -----------------------------------------------------
def system_call(command):
  """
  Makes a system call
  """
  info("  External command executed: "+command)
  #subprocess.call(command.split())
  os.system(command)

# -----------------------------------------------------
def get_time():
  """
  Gets the time in human-readable format
  """
  return time.asctime(time.gmtime())

# -----------------------------------------------------
def info(text):
  """
  Prints an info into the log-file
  """
  msg = get_time() + ': '+text
  
  log_file = cp.get('paths','main')+'/llmonitor.log'
  logfile = open(log_file,'a')
  logfile.write(msg+'\n')
  logfile.close()


# -----------------------------------------------------
def send_mail(subject, msg):
  """
  Function to send an email to a certain adress
  """

  # Adjust messages and subjects automatically
  message = 'Automatic notification from pylal_exttrig_llmonitor at time '+\
            get_time()+'\n\n'+subject+'\n'+msg
  subject = 'S6-ExtOnline: '+subject
    
  # open file for detailed output message
  tmp_file = '.llmonitor.email'
  f = file(tmp_file,'w')
  f.write(message)
  f.close()

  # call the command
  email_adresses = cp.get('notifications','email').replace(',',' ').split()
  for address in email_adresses:
    command = "mail -s '%s' %s < %s" % (subject, address, tmp_file)
    system_call(command)
  

# -----------------------------------------------------
def notify(dag, message):
  """
  Makes an email notification to all recipients listed
  in the config file
  """

  # construct the subject of the email
  subject = 'Status changed for DAG GRB%s: %s' %\
       (dag['name'], message)

  # construct the message for the email
  email_msg = 'Automatic notification from pylal_otex at time %s\n\n'%\
              get_time()
  email_msg += subject+'\n'
  email_msg += 'Stage currently executed: %s\n'%dag['stage']
  email_msg += 'The DAG is located at : %s\n'% dag['path']

  # send the email to all recipients
  send_mail(subject, message)

  # and note it in the log-file
  info("  Email notification sent with the following content: "+\
       email_msg.replace('\n','\n    '))

# --------------------------------------
def get_dag_part(ini_file):
  """
  Gets the dag-name from the ini file.
  This might be non-robust, therefore it is
  coded as a complete function.
  """

  dag_part = ini_file.split('.')[0]
  return dag_part

# --------------------------------------
def get_segment_lists():
  """
  Function to download the latest segment lists
  """
  


# --------------------------------------
def get_empty_exttrig_row():
  """
  Create an empty exttrig row 
  """
  row = lsctables.ExtTriggersTable()

  row.process_id = None
  row.det_alts = None
  row.det_band = None
  row.det_fluence = None
  row.det_fluence_int = None
  row.det_name = None
  row.det_peak = None
  row.det_peak_int = None
  row.det_snr = ''
  row.email_time = 0
  row.event_dec = 0.0
  row.event_dec_err = 0.0
  row.event_epoch = ''
  row.event_err_type = ''
  row.event_ra = 0.0
  row.event_ra_err = 0.0
  row.start_time = 0
  row.start_time_ns = 0
  row.event_type = ''
  row.event_z = 0.0
  row.event_z_err = 0.0
  row.notice_comments = ''
  row.notice_id = ''
  row.notice_sequence = ''
  row.notice_time = 0
  row.notice_type = ''
  row.notice_url = ''
  row.obs_fov_dec = 0.0
  row.obs_fov_dec_width = 0.0
  row.obs_fov_ra = 0.0
  row.obs_fov_ra_width = 0.0
  row.obs_loc_ele = 0.0
  row.obs_loc_lat = 0.0
  row.obs_loc_long = 0.0
  row.ligo_fave_lho = 0.0
  row.ligo_fave_llo = 0.0
  row.ligo_delay = 0.0
  row.event_number_gcn = 0
  row.event_number_grb = ''
  row.event_status = 0
  return row

# --------------------------------------
def generate_summary(monitor_file):

  def add(table, text):
    #print text
    return table + '<td>' +str(text)+'</td>' 

  # Read the list of all processed GRBs
  monitor_list = pickle.load(file(monitor_file))

  # Bring them in timely order
  time_unsort = [grb['triggertime'] for grb in monitor_list]
  index = np.argsort(time_unsort)

  table = total_summary_prefix % get_time()

  # Loop over all GRBs in reverse order
  for i in index[::-1]:

    item = monitor_list[i]
    name = item['name']

    if (i % 2):
      table += '<tr style="background-color: rgb(255, 200, 200);">'
    else:
      table += '<tr style="background-color: rgb(255, 200, 200);">'

    table = add(table, i+1)
    table = add(table, name) 
    table = add(table, item['status'])
    table = add(table, item['triggertime'])
    asctime = time.asctime(time.gmtime(item['triggertime']+offset_gps_to_linux))
    table = add(table, asctime)
    if item.has_key('redshift'):
      table = add(table, '<a href="http://gcn.gsfc.nasa.gov/gcn3/%d.gcn3">%.2f</a>' % (item['redshift_ref'], item['redshift']))
    else:
      table = add(table, '&mdash')
    if item.has_key('duration'):
      table = add(table, '<a href="http://gcn.gsfc.nasa.gov/gcn3/%d.gcn3">%.2f</a>' % (item['duration_ref'], item['duration']))
    else:
      table = add(table, '&mdash')
    table = add(table, '%.2f' % item['right_ascension'])
    table = add(table, '%.2f' % item['declination'])
    for ifo in ifo_list:
      if ifo in item['ifos']:
        table = add(table, '<b>%.2f</b>'%item['qvalues'][ifo])
      else:
        table = add(table, '%.2f'%item['qvalues'][ifo])


    if item['status'] == 'Finished':
      table = add(table, 'RESULTS')
      htmlfile = cp.get('paths','publishing_url')+'GRB%s/pylal_exttrig_llsummary_%s.html' % (name, name)
      table = add(table, '<a href="%s">link</a>'%htmlfile)
    else:
      table = add(table, '&mdash')
      table = add(table, '&mdash')

    table +='</tr>'

    # Processing only the one having no data of having finished
    #if item['stage']=='NoData' or item['stage']=='Finished'
    #table += ""
    
    filename = cp.get('paths','publishing_path')+'/total_summary.html'
    f = open(filename,'w')
    f.write(table)
    f.close()

# -----------------------------------------------------
class ExttrigDag(object):
  """
  Class holding all the infos for a new analysis DAG.
  This has the advantage that the DAG creating stuff
  can be more easily split up intop a couple of shorter
  functions.
  """

  # -----------------------------------------------------
  def __init__(self, grb_name=None, grb_ra=None, grb_de=None, grb_time=None):
    self.name = grb_name # just the number, i.e. 091023C
    self.ra = float(grb_ra)
    self.de = float(grb_de)
    self.time = int(grb_time)

    self.qvalues = {}

  # -----------------------------------------------------
  def set_paths(self, input_dir=None, glue_dir=None, pylal_dir = None,\
                lalapps_dir=None, main_dir=None,\
                condor_log_path = None, log_file=None):
    self.input_dir = input_dir
    self.glue_dir = glue_dir
    self.pylal_dir = pylal_dir
    self.lalapps_dir = lalapps_dir
    self.main_dir = main_dir
    self.condor_log_path = condor_log_path
    self.log_file = log_file
     
    self.analysis_dir = self.main_dir+'/GRB'+self.name

  # -----------------------------------------------------
  def set_ini_file(self, ini_file):
    self.ini_file = ini_file
    self.dag_file = get_dag_part(self.ini_file)+'_uberdag.dag'
    self.out_file = self.dag_file+'.dagman.out'

  # -----------------------------------------------------
  def set_monitor_file(self, monitor_file):
    self.monitor_file = monitor_file

  # -----------------------------------------------------
  def set_addresses(self, addresses):
    self.addresses = addresses

  # -----------------------------------------------------
  def set_seg_info(self, offsource_segment = None, ifos = None):
    self.offsource_segment = offsource_segment
    self.ifos = ifos

  def make_links(self, sourcedir, destdir, list_exec):
    """
    Using a list of executables names, make symbolic links rather
    than copying the files itself. That way the origin, i.e. the 
    version of the executables can be verified, at least somewhat...
    """
    
    cmd = 'cd %s;' % destdir
    for execfile in list_exec:
      cmd += 'ln -s %s/%s .;' % (sourcedir, execfile)
    system_call(cmd)


  # -----------------------------------------------------
  def create_exttrig_xml_file(self):
    """
    Creates an exttrig xml file with the basic informations
    of the GRB in it.
    Data given in strings or as float
    """
  
    #
    # First, create a dummy XML file with the informations in it
    #
    
    # prepare a new XML document
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    tbl = lsctables.New(lsctables.ExtTriggersTable)
    xmldoc.childNodes[-1].appendChild(tbl)
    
    # set the values we need
    row = get_empty_exttrig_row()
    row.event_ra = float(self.ra)
    row.event_dec = float(self.de)
    row.start_time = int(self.time)
    row.event_number_gcn = 9999
    row.event_number_grb = self.name
    
    # insert into the table
    tbl.extend([row])
    
    # write out the trigger file
    self.trigger_file = 'GRB%s.xml' % self.name
    utils.write_filename(xmldoc, self.trigger_file)

  # -----------------------------------------------------
  def run_datafind(self):
    """
    Run the datafind command to find the data
    """ 
    # create the datafind directory 
    cache_dir = "%s/GRB%s/datafind/cache" % (self.analysis_dir, self.name)
    cmd = 'mkdir -p '+cache_dir
    system_call(cmd)

    executable = self.lalapps_dir+'/bin/lalapps_online_datafind'
    #print self.offsource_segment
    #print executable
    starttime = self.offsource_segment[0]
    endtime = self.offsource_segment[1]

    # and run the datafind command for each IFO, putting
    # the cache files directly into them
    for ifo in ifo_list:
      if ifo=='V1':
        cache_filename = "V-V1_DMT_HREC-%9d-%9d.cache" % (starttime, endtime)
      else:
        cache_filename = "%s-DMT_C00_L2-%9d-%9d.cache" %(ifo[0].upper(), starttime, endtime)

      cmd = "%s --ifo %s --gps-start-time %d --gps-end-time %d --output %s/%s" % \
            (executable, ifo, starttime, endtime, cache_dir, cache_filename)
      system_call(cmd)


  # -----------------------------------------------------
  def prepare_analysis_directory(self, skip_rundag = False):
    #
    # Now create directories and copy a bunch of files
    #
    
    # Create the main directory
    if False and os.path.exists(self.analysis_dir):
      raise IOError, "The directory %s already exists. Please (re)move"\
            " this directory or choose another name for the "\
            "analysis directory" % self.analysis_dir
    cmd = 'mkdir %s' % self.analysis_dir
    system_call(cmd)

    # copy some basic files into this directory
    cmd = 'cp %s/*.ini %s/' % \
         (self.input_dir, self.analysis_dir)
    system_call(cmd)
    cmd = 'cp %s %s/' % (self.trigger_file, self.analysis_dir)
    system_call(cmd)

    # replace some arguments in the ini-file
    config_file = '%s/%s' % (self.analysis_dir, self.ini_file)
    pc = ConfigParser.ConfigParser()
    pc.read(config_file)
    pc.set("pipeline", "user-tag", "GRB%s" %self.name)

    # write out new ini-file
    cp_file = open(config_file, 'w')
    pc.write(cp_file)
    cp_file.close()


    # copy executables
    cmd = 'cd %s/bin; cp lalapps_coherent_inspiral lalapps_coherentbank \
      lalapps_coire lalapps_frjoin lalapps_inca lalapps_inspinj lalapps_inspiral \
      lalapps_inspiral_hipe lalapps_sire lalapps_thinca lalapps_tmpltbank \
      lalapps_trigbank lalapps_plot_hipe lalapps_trigger_hipe %s' %\
    (self.lalapps_dir, self.analysis_dir)
    system_call(cmd)

    #self.make_links(self, self.glue_dir+'/bin', self.analysis_dir, ['ligo_data_find','ligolw_add'])
    cmd = 'cd %s/bin; cp ligo_data_find ligolw_add %s' % \
      (self.glue_dir, self.analysis_dir)
    system_call(cmd)

    cmd = 'cd %s/bin; cp pylal_relic pylal_exttrig_llsummary %s' % \
      (self.pylal_dir, self.analysis_dir)
    system_call(cmd)


    # copy the segment files
    cmd = 'cp %s/*-science.txt %s' % (self.main_dir, self.analysis_dir)
    system_call(cmd)

    #
    # make the call to trigger_hipe
    #
    cmd = 'cd %s;' % self.analysis_dir
    cmd += template_trigger_hipe % \
           (self.trigger_file, self.name, self.ini_file, self.condor_log_path)
    system_call(cmd)

    # TODO HERE
    self.run_datafind()

    # Prepare the postprocesing directory at this stage
    path = "%s/GRB%s/postprocessing" % (self.analysis_dir, self.name)
    system_call('mkdir -p %s/logs'%path)

    cmd = 'cp %s/post* %s' % (self.input_dir, path)
    system_call(cmd)


  # -----------------------------------------------------
  def start_analysis(self):

    # create the call to start the DAG
    cmd = 'cd %s;' % self.analysis_dir
    cmd += 'export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR=FALSE;'
    cmd += 'condor_submit_dag %s' % self.dag_file
    if skip_rundag:
      info("Execution of condor_submit skipped. Command would have been: "+cmd)
    else:
      system_call(cmd)

  # -----------------------------------------------------
  def check_analysis_directory(self):
    """
    Check if the dagman.out file does exist
    after some while
    """

    out_file = self.analysis_dir+'/'+self.out_file

    success = False
    for i in range(10):
      time.sleep(10)

      # try to open the file
      try:
        f = open(out_file)
        f.close()
        # if no errors: great, DAG is running...
        success = True
        break
      except:
        success = False

    # check the status at the end
    if not success:
      
      # send an email about this problem
      subject = 'Problems starting condor DAG'     
      email_msg = 'The condor DAG %s was not started' % self.dag_file
      email_msg += 'The DAG is located at : %s\n'% self.analysis_dir
      send_mail(subject, email_msg)  
      
      return -1 
    else:

      return 1
      
  # -----------------------------------------------------
  def update_database(self, has_data):
    """
    Update the database, put in the informations
    for this GRB:
    """

    #
    # add this DAG to the monitor database
    #

    dag_dict = {'name':self.name, 'path':self.analysis_dir, \
                 'dag-name':self.dag_file,\
                 'ifos':'', 'condorlogpath':'',\
                 'triggertime':self.time,\
                 'right_ascension': self.ra,'declination':self.de,
                 'qvalues':self.qvalues}

    if has_data ==0:
      # Not enough data has been found. The data is stored,
      # but no analysis is started
      dag_dict['starttime'] = None
      dag_dict['endtime'] = None
      dag_dict['status'] = 'NoData'
      dag_dict['stage'] = 'NoData'
      dag_dict['ifos'] = ''
      dag_dict['condorlogpath'] = ''
    else:
      # Enough data is stored, the DAG has been prepared
      dag_dict['starttime'] = self.offsource_segment[0]
      dag_dict['endtime'] = self.offsource_segment[1]
      # The DAG is either running or has not started.
      # The status of the DAG will be checked then in the next round
      if has_data==1:
        dag_dict['status'] = 'Running'
      if has_data ==-1:
        dag_dict['status'] = 'Aborted'
      dag_dict['stage'] = 'Inspiral' 
      dag_dict['ifos'] = self.ifos
      dag_dict['condorlogpath'] = self.condor_log_path

    # update the monitor file
    monitor_list = pickle.load(file(self.monitor_file))
    monitor_list.append(dag_dict)
    pickle.dump(monitor_list, file(self.monitor_file,'w'))

    if has_data:
      info('The GRB %s has been added to the database'%self.name)
    else:
      info('The GRB %s has been added to the database with NoData.'%self.name)


  # -----------------------------------------------------
  def get_minimum_scienceseg_length(self):
    """
    Calculate the minimum science segment that
    can be used with the data given in the actual ini-file.
    The procedure below is from trigger_hipe, after
    having replaced 'cp by 'pc'
    """
    
    # get the name of the ini-file to be used
    ini_file = cp.get('paths','cvs') + '/'+cp.get('data','ini_file')
 
    # the following is just a copy-and-paste from trigger_hipe
    # having replaced 'cp' by 'pc'
    pc = ConfigParser.ConfigParser()
    pc.read(ini_file)
    paddata = int(pc.get('data', 'pad-data'))
    n = int(pc.get('data', 'segment-length'))
    s = int(pc.get('data', 'number-of-segments'))
    r = int(pc.get('data', 'sample-rate'))
    o = int(pc.get('inspiral', 'segment-overlap'))
    length = ( n * s - ( s - 1 ) * o ) / r
    overlap = o / r
    minsciseg = length + 2 * paddata
    
    return minsciseg

  # --------------------------------------
  def update_segment_lists(self, timeoffset):
    """
    Function to download the latest segment lists.
    @param timeoffset: The offset in time for downloading those segments
    """

    seg_names = ['H1:DMT-SCIENCE:1','L1:DMT-SCIENCE:1','V1:ITF_SCIENCEMODE']
    ifo_list = ['H1','L1','V1']
    starttime = self.time-timeoffset
    endtime = self.time+timeoffset

    for ifo, seg in zip(ifo_list, seg_names):

      segxmlfile = "%s/segments%s.xml" % (self.main_dir, ifo)
      segtxtfile = "%s/%s-science.txt" % (self.main_dir, ifo)

      cmd = "%s/bin/ligolw_segment_query --database --query-segments --include-segments '%s' --gps-start-time %d --gps-end-time %d "\
                  "> %s" %\
                  (self.glue_dir, seg, starttime, endtime, segxmlfile)
      system_call(cmd)

      # 'convert' the data from the xml format to a useable format...
      # TODO: change the other places to accept the xml format
      doc = utils.load_filename(segxmlfile)
      segs = table.get_table(doc, "segment")
      seglist = segments.segmentlist(segments.segment(s.start_time, s.end_time) for s in segs)
      segmentsUtils.tosegwizard(file(segtxtfile, 'w'), seglist, header = True)
  
  # -----------------------------------------------------
  def get_segment_info(self,plot_segments_file = None):

    minsciseg = self.get_minimum_scienceseg_length()
    segdict = segments.segmentlistdict()

    # get the segment dicts
    for ifo in ifo_list:
      ifo_segfile = '%s/%s-science.txt' % (self.main_dir, ifo)
      if ifo_segfile is not None:
        tmplist = segmentsUtils.fromsegwizard(open(ifo_segfile))
        segdict[ifo] = segments.segmentlist([s for s in tmplist \
                                             if abs(s) > minsciseg])
    ifolist = segdict.keys()
    ifolist.sort()

    # create the onsource segment
    onsource_left = int(cp.get('data','onsource_left'))
    onsource_right = int(cp.get('data','onsource_right'))
    trigger = int(self.time)
    onSourceSegment = segments.segment(trigger - onsource_left,
                                       trigger + onsource_right)

    # convert string in integer
    padding_time = int(cp.get('data','padding_time'))
    num_trials = int(cp.get('data','num_trials'))
    symmetric = False
    offSourceSegment, grb_ifolist = micos(segdict, onSourceSegment,
                                          padding_time = padding_time, \
                                          max_trials = num_trials,
                                          min_trials = num_trials, \
                                          symmetric = symmetric)


    grb_ifolist.sort()
    ifo_times = "".join(grb_ifolist)

    # make a plot of the segments if required
    if plot_segments_file:

      plot_offset = 1000

      length_off_source = num_trials*(abs(onSourceSegment))
      plot_offSourceSegment = segments.segment(onSourceSegment[0] - length_off_source,
                                          onSourceSegment[1] + length_off_source)

      effective_window = segments.segmentlist([plot_offSourceSegment]).\
                         protract(plot_offset)
      effective_segdict = segdict.map(lambda sl: sl & effective_window)
      plot = PlotSegmentsPlot(trigger)
      plot.add_contents(effective_segdict)
      plot.set_window(offSourceSegment, plot_offset)
      plot.highlight_segment(onSourceSegment)
      plot.finalize()
      plot.ax.set_title('Segments for GRB '+self.name)
      plot.savefig(plot_segments_file)
      plot.close()

    # return the main values from this function    
    return offSourceSegment, grb_ifolist

  # -----------------------------------------------------
  def calculate_optimality(self):

    # Calculate the antenna factors (for informational purpose only)
    for ifo in ifo_list:
      _, _, _, q = antenna.response(self.time, self.ra, self.de, \
                                    0.0, 0.0, 'degree', ifo)
      self.qvalues[ifo]=q

