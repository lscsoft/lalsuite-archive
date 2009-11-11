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
from datetime import datetime

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
from pylal.xlal import date

# the config parser to be used in some of the functions
cp = None

template_trigger_hipe = "./lalapps_trigger_hipe"\
  " --h1-segments H1-science_grb%s.txt" \
  " --l1-segments L1-science_grb%s.txt" \
  " --v1-segments V1-science_grb%s.txt" \
  " --list %s" \
  " --grb %s --onsource-left 5" \
  " --onsource-right 1 --config-file %s" \
  " --injection-config injectionsWI.ini --log-path %s" \
  " --number-buffer-left 8 --number-buffer-right 8" \
  " --num-trials 340 --padding-time 72 --verbose" \
  " --skip-datafind --skip-dataquality --user-tag onoff"

template_trigger_hipe_inj = "./lalapps_trigger_hipe"\
  " --h1-segments H1-science_grb%s.txt" \
  " --l1-segments L1-science_grb%s.txt" \
  " --v1-segments V1-science_grb%s.txt" \
  " --list %s" \
  " --grb %s --onsource-left 5" \
  " --onsource-right 1 --config-file %s" \
  " --injection-config %s --log-path %s" \
  " --number-buffer-left 8 --number-buffer-right 8" \
  " --num-trials 340 --padding-time 72 --verbose" \
  " --skip-datafind --skip-dataquality --skip-onsource --skip-offsource"\
  " --user-tag inj --overwrite-dir"


ifo_list = ['H1','L1','V1']

offset_gps_to_linux = 315964800 # see http://www.epochconverter.com/ for 6 Jan 1980 00:00:00 GPS = 000000000

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
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Status OO</td>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Status INJ</td>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">GPS<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Date<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">redshift<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">duration [s]<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">RA<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">DEC<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">H1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">L1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">V1<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Sanity<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Result<br>
  <td style="vertical-align: top; font-weight: bold; font-style: italic; color: rgb(51, 51, 255); background-color: rgb(255, 153, 0);">Box<br>
"""

# -----------------------------------------------------
def system_call(item, command, divert_output_to_log = True):
  """
  Makes a system call.
  @params item: a text specifying the content of the text
         (e.g. number of the GRB the message is associated with)
         (see also 'info')
  @params command: the command to be executed on the bash
  @params divert_output_to_log: If this flag is set to True the output of the 
                                given command is automatically put into the log-file.
                                If the output of some command itself is further used,
                                like science segments, this flag must be set 
                                to False, so that the output is diverted where it should go.
  """
  l = logfile_name()

  # put the command used into the log file
  info(item, ">>> "+command)
  
  # and the output (and error) of the command as well
  #subprocess.call(command+' >%s 2>%s '%(l,l), shell=True)
  if divert_output_to_log:
    os.system(command +' >>%s 2>>%s '%(l,l))
  else:
    os.system(command +' 2>>%s '%l)

# -----------------------------------------------------
def get_time():
  """
  Returns the current time in human-readable format
  """
  return time.asctime(time.gmtime())

# -----------------------------------------------------
def get_gps_from_asc(date_string, time_string):
  """
  Computes the correct GPS time from the date and time
  as given in text strings.
  @param date_string: date in a string format, i.e. 090717
  @param time_string: time in a string format, i.e. 19:10:34
  """

  # convert the date and times (as read from the trigger file)
  # into tuples
  a = time.strptime(date_string, "%y%m%d")
  time_list = time_string.split('.')
  b = time.strptime(time_list[0], "%H:%M:%S")
  if len(time_list)==2:
    nsecs = time_list[1]
    nsecs += (9-len(nsecs))*'0'
    nano_seconds = int(nsecs) 
  else:
    nano_seconds = 0

  # populate a datetime tuple
  tm = datetime(a[0], a[1], a[2], b[3], b[4], b[5]).timetuple()
  # and parse it, the last three entries populated as well,
  # to the wrapped XLALUTCToGPS function.
  gpstime = date.XLALUTCToGPS(tm)

  return int(gpstime)


# -----------------------------------------------------
def logfile_name():
  """
  Returns the file of the logfile; used in 'info' and 'system_call'
  """
  return cp.get('paths','main')+'/llmonitor.log'

# -----------------------------------------------------
def info(item, text):
  """
  Prints an info into the log-file.
  @item: a text specifying the content of the text 
         (e.g. number of the GRB the message is associated with)
  @text: the text to be logged
  """
  msg = get_time() + ' ('+item+'): '+text
  
  log_file = logfile_name()
  logfile = open(log_file,'a')
  logfile.write(msg+'\n')
  logfile.close()

  print text

# -----------------------------------------------------
def send_mail(subject, msg, email_adresses = None):
  """
  Function to send an email to a certain adress
  @param subject: Subject line of the email
  @param msg: Message of the email
  @param email_adresses: list of email adresses to which the mail is sent
  """

  # Adjust messages and subjects automatically
  message = 'Automatic notification from pylal_exttrig_llmonitor at time '+\
            get_time()+'\n\n'+subject+'\n'+msg
  subject = cp.get('notifications','head') + ': '+subject
    
  # open file for detailed output message
  tmp_file = '.llmonitor.email'
  f = file(tmp_file,'w')
  f.write(message)
  f.close()

  # select the recipients
  if not email_adresses:
    email_adresses = cp.get('notifications','email').replace(',',' ').split()

  # send the message to all recipients
  for address in email_adresses:
    command = "mail -s '%s' %s < %s" % (subject, address, tmp_file)
    system_call('email',command)
 
  #info('email','Email content: '+msg) 

# -----------------------------------------------------
def notify(grb, dag, message):
  """
  Makes an email notification to all recipients listed
  in the config file.
  @param grb: grb dictionary for obtaining some informations
  @param message: the message of the notification
  """

  # construct the subject of the email
  subject = 'Status changed for DAG GRB%s: %s' %\
       (grb.name, message)

  # construct the message for the email
  email_msg = 'Automatic notification from pylal_exttrig_llutils at time %s\n\n'%\
              get_time()
  email_msg += subject+'\n'
  email_msg += 'The analysis dir is %s\n' % grb.analysis_dir
  email_msg += ' and the dagfils is %s\n' % dag.get_outname()

  # send the email to all recipients
  send_mail(subject, message)

  # and note it in the log-file
  info("email","  Email notification sent with the following content: "+\
       email_msg.replace('\n','\n    '))


# --------------------------------------
def get_lockname():
  """
  Returns the name of the lock file
  """
  return cp.get('paths','main')+'/.llmonitor.lock'

# --------------------------------------
def check_lock():
  """
  Checks if the lock file exists
  """
  return os.path.exists(get_lockname())

# --------------------------------------
def set_lock():
   """
   Sets the lock file
   """
   f = open(get_lockname(),'w')
   f.write(get_time())
   f.close()

# --------------------------------------
def del_lock():
   """
   Removes the lock file
   """
   os.remove(get_lockname())


# --------------------------------------
def get_dag_part(ini_file):
  """
  Gets the dag-name from the ini file.
  This might be non-robust, therefore it is
  coded as a complete function which can be changes easily.
  @param ini_file: the name of the ini-file
  @param return: the common part of any dag name
  """
  dag_part = ini_file.split('.')[0]
  return dag_part

# --------------------------------------
def get_empty_exttrig_row():
  """
  Returns an empty exttrig row 
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
def get_monitor_filename():
  """
  Returns the name of the monitor pickle filename
  """
  return cp.get('paths','main')+'/llmonitor.pickle'

# --------------------------------------
def read_monitor_list():
  """
  Opens the monitor pickle file (usually llmonitor.pickle)
  and return its contents
  """

  monitor_file = get_monitor_filename()
  try:
    monitor_list = pickle.load(file(monitor_file))
  except IOError:
    # create an empty file if it does not exist
    monitor_list = []
    pickle.dump(monitor_list, file(monitor_file,'w'))
  return monitor_list

# --------------------------------------
def write_monitor_list(monitor_list):
  """
  Writes the monitor list to file
  @param monitor_list: list to be written to file
  """
  monitor_file = get_monitor_filename()
  pickle.dump(monitor_list, file(monitor_file,'w'))

# --------------------------------------
def read_grb_from_list(grb_name):
  """
  Returns the object associated with the given GRB.
  @params grb_name: name of the GRB without the leading 'GRB'
  """
  grb_list = read_monitor_list()
  for grb in grb_list:
    if grb.name==grb_name:
      return grb
  return None

# --------------------------------------
def copy_exttrig_nofications():
  """
  Copying all relevant files to the working directory,
  usually from CIT from Isabels pwd
  """
  alert_loc = cp.get('alerts','alert_location')
  main_loc = cp.get('paths','main')
  cmd = 'scp %s %s >> ~/cp.log 2>&1' % (alert_loc, main_loc)
  system_call('monitor', cmd)

# --------------------------------------
def update_durations(monitor_list):
  """
  Reads the local copy of the parsed circular and 
  updated any duration information in the monitor_list structure
  @params monitor_list: list of all GRBs and DAGs
  """
  # name of the circular file
  circular_file = cp.get('paths','main')+'/'+cp.get('alerts','circular_file')

  # Read all durations from the circular file
  dict_duration = {}
  for line in file(circular_file):
    parts = line.split()
    grb_name = parts[2]
    duration = float(parts[12])

    # store the duration. A value of zero means it is unknown
    if duration>0:
      dict_duration[grb_name]=duration

  # Read the list of all processed GRBs
  for grb in monitor_list:
    
    # update the duration information when available
    if grb.name in dict_duration:
      grb.duration = dict_duration[grb.name]


# --------------------------------------
def obtain_results(grb):
  """
  Obtain the result, i.e the smallest p(c|0)
  @param grb: the grb stucture with all the infos in it
  """

  name = grb['name']
  path_to_result = '%s/GRB%s/postprocessing/output/llsummary_GRB%s.pickle' %\
     (grb['path'], name, name)
  data = pickle.load(file(path_to_result))

  min_prob = 2.0
  for coinc in data:
    if 'prob' in coinc:
      p = coinc['prob']
      if p<min_prob:
        min_prob = p

  return min_prob

# --------------------------------------
def generate_summary(publish_path, publish_url):
  """
  Generating summary page, with all sanity and/or openbox results
  properly linked.
  @param publish_path: Main path to where to copy the results and files
  @param publish_url: The url identifier of the same path
  """

  def add(table, text):
    return table + '<td>' +str(text)+'</td>' 

  def add_linked_value(table, value, ref):
    if value>0:
      if ref>0:
        table = add(table, '<a href="http://gcn.gsfc.nasa.gov/gcn3/%d.gcn3">%.2f</a>' % (ref, value))
      else:
        table = add(table, '%.2f' % value)
    else:
      table = add(table, '&mdash')
    return table


  # Read the list of all processed GRBs
  monitor_list = read_monitor_list()

  # Bring them in timely order
  time_unsort = [grb.time for grb in monitor_list]
  index = np.argsort(time_unsort)
  num_grb = len(time_unsort)

  table = total_summary_prefix % get_time()

  # Loop over all GRBs in reverse order
  for number, i in enumerate(index[::-1]):

    grb = monitor_list[i]

    # make the table background coloring
    if grb.duration and grb.duration<float(cp.get('data','max-duration')):
      if grb.has_data:
        coldef = '255, 200, 200'
      else:
        coldef = '180, 180, 180'
    else:
      if grb.has_data:
        coldef = '153, 255, 255'
      else:
        coldef = '150, 150, 150'

    table += '<tr style="background-color: rgb(%s);">' % coldef
    
    if grb.has_data:
      status_onoff = grb.dag['onoff'].get_status()
      status_inj = grb.dag['inj'].get_status()
    else:
      status_onoff = status_inj = 0
    table = add(table, num_grb- number)
    table = add(table, '<a href="http://grblog.org/grblog.php?view=burst&GRB=%s">%s</a>'%(grb.name, grb.name)) 
    status_msg = grb.get_html_status()
    table = add(table, status_msg['onoff'])
    table = add(table, status_msg['inj'])
    table = add(table, grb.time)
    asctime = time.asctime(time.gmtime(grb.time+offset_gps_to_linux))
    table = add(table, asctime)
    table = add_linked_value(table, grb.redshift, None )
    table = add_linked_value(table, grb.duration, None)
    table = add(table, '%.2f' % grb.ra)
    table = add(table, '%.2f' % grb.de)
    for ifo in ifo_list:
      segplot_link = 'GRB%s/plot_segments_grb%s.png'%(grb.name, grb.name)
      
      if ifo in grb.ifos:
        txt = '<b>%.2f</b>'%grb.qvalues[ifo]
      else:
        txt = '%.2f'%grb.qvalues[ifo]
      table = add(table, '<a href="%s">%s</a>'%(segplot_link, txt))
    
    if status_onoff==5:
     
      #Add link to sanity pages
      htmlfile = publish_url+'/GRB%s/pylal_exttrig_llsummary_%s-sanity.html' % (grb.name, grb.name)
      htmlfile_inj = publish_url+'/GRB%s/pylal_exttrig_llsummary_%s-sanity_inj.html' % (grb.name, grb.name)
      if status_inj==5:
        table = add(table, '<a href="%s">onoff</a> <a href="%s">inj</a> '%(htmlfile, htmlfile_inj))
      else: 
        table = add(table, '<a href="%s">onoff</a> &mdash '%htmlfile)

      # Add link to box
      if grb.openbox:
        # add result    
        result = obtain_results(grb)
        if result<2:
          table = add(table, '%.2f'%result)
        else:        
          table = add(table, 'no cand.')

        # and link to the openbox details
        htmlfile = publish_url+'/GRB%s/pylal_exttrig_llsummary_%s.html' % (grb.name, grb.name)
        htmlfile_inj = publish_url+'/GRB%s/pylal_exttrig_llsummary_%s_inj.html' % (grb.name, grb.name)
        if status_inj==5:
          table = add(table, '<a href="%s">onoff</a> <a href="%s">lik</a> '%(htmlfile, htmlfile_inj))
        else:
          table = add(table, '<a href="%s">onoff</a> &mdash '%htmlfile)

      else:
        # box closed otherwise
        table = add(table,'box closed')
        table = add(table,'box closed')
    else:
      # no data available if analysis not finished
      table = add(table, '&mdash')
      table = add(table, '&mdash')
      table = add(table, '&mdash')

    table +='</tr>'

    # write out the complete html file
    filename = publish_path +'/total_summary.html'
    f = open(filename,'w')
    f.write(table)
    f.close()

# -----------------------------------------------------
# -----------------------------------------------------
class AnalysisDag(object):
  """
  Class to hold and handle an analysis DAG and all
  related information.
  """
  
  # -----------------------------------------------------
  def __init__(self, name, type, analysis_dir):
    """
    Initializing this class with all the needed information
    @param name: name of the GRB
    @param type: what dag is this? onoff/inj
    #@param stage: stage of the dag, like uberdag or ligolwdag
    #@param inifile: inifile for this DAG
    #@param injfile: injection file for this DAG
    @param analysis_dir: path to the analysis directory
    """

    # store the input data
    self.name = name
    self.type = type
    self.analysis_dir = analysis_dir

    self.dagname = None

    self.status = 0
    self.status_dict = {1:'inspiral',2:'ligolw',3:'postproc'}

  # --------------------------------------
  def set_dagname(self, name):
    """
    Sets the current name of the DAG
    @param name: name of the .dag file
    """
    self.dagname = name

  # --------------------------------------
  def get_basename(self):
    """
    returns basename without any ending.
    """
    return self.analysis_dir+'/'+self.dagname[:-4]
  
  # --------------------------------------
  def get_outname(self):
    """
    Returns the outname of this DAG
    """
    return self.dagname+'.dagman.out'

  # --------------------------------------
  def get_dagname(self):
    """
    Returns the name of the DAG file
    """
    return self.dagname
 
  # --------------------------------------
  def get_shname(self):
    """
    Returns the name of the sh file
    """
    return self.get_basename()+'.sh'

  # --------------------------------------
  def start(self):
    """
    Start this DAG
    """

    # create the call to start the DAG
    cmd = 'cd %s;' % self.analysis_dir
    cmd += 'export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR=FALSE;'
    cmd += 'condor_submit_dag %s' % self.get_dagname()
    #print "WARNING: DAG NOT STARTED FOR TESTING"
    system_call(self.name, cmd)

    # change the status
    self.status = 1

  # --------------------------------------
  def get_status(self):
    return self.status

  # --------------------------------------
  def set_status(self, new_status):
    if new_status<=0:
      raise ValueError, "The DAG Status variable can only be set to positive values"

    self.status = new_status

  # --------------------------------------
  def get_stage_name(self):
    """
    Returns the name of the stage or error of the current DAG.
    """
    status_dict = {1:'inspiral',2:'ligolw',3:'postproc'}
    # check the status number and choose the text
    text = ''
    if self.status==0:
      text = 'Not started'
    elif self.status==5:
      text = 'Complete'
    elif self.status==-6:
      text = 'DAGFILE ERROR'
    else:
      text = status_dict[abs(self.status)]
      if self.status<0:
        text += "ERROR"
 
    return text

  # --------------------------------------
  def check_status(self, grb, dag):
    """
    Updating the status for this DAG,
    and return the fstat value
    """

    # try to open the dagman.out file
    fstat = 0
    try:
      # read the last line only
      line = file(self.get_outname()).readlines()[-1]
      if "EXITING WITH STATUS 1" in line:
        fstat = -1
      if "EXITING WITH STATUS 0" in line:
        fstat = 1
    except IOError:
      fstat = -2

    # change the status if the DAG was running before
    if self.status>0:
     if fstat<0:
       # set the status to error
       self.status = -self.status

       if fstat == -1:
         notify(grb, dag, 'DAG exited on error')
       elif fstat==-2:
         notify(grb, dag, 'DAG file vanished!?')

    return fstat

  
# -----------------------------------------------------
# -----------------------------------------------------
class GRB(object):
  """
  Class holding all the infos for a GRB and for setting up
  a new analysis DAG.
  """

  # -----------------------------------------------------
  def __init__(self, grb_name=None, grb_ra=None, grb_de=None, grb_time=None):
    self.name = grb_name # just the number, i.e. 091023C
    self.ra = float(grb_ra)
    self.de = float(grb_de)
    self.time = int(grb_time)

    self.ifos = ''
    self.duration = None
    self.redshift = None
    self.starttime = None
    self.endtime = None
    self.openbox = False
    
    # prepare the DAG instances
    self.dag = {'onoff':None, 'inj':None}

    # prepare variables for later use
    self.qvalues = {}
    self.offsource_segment = None
    self.ifolist  = []

    self.use_offline_data = False
    self.type_online = {'H1':'H1_DMT_C00_L2', 'L1':'L1_DMT_C00_L2', 'V1':'V1_DMT_HREC'}
    self.type_offline = {'H1':'H1_LDAS_C00_L2','L1':'L1_LDAS_C00_L2','V1':'HrecOnline'}

  # -----------------------------------------------------
  def set_paths(self, input_dir=None, glue_dir=None, pylal_dir = None,\
                lalapps_dir=None, main_dir=None,\
                ini_file = None, inj_file = None,\
                config_file = None, \
                condor_log_path = None, log_file=None):
    self.input_dir = input_dir
    self.glue_dir = glue_dir
    self.pylal_dir = pylal_dir
    self.lalapps_dir = lalapps_dir
    self.main_dir = main_dir
    self.inifile = ini_file
    self.injfile = inj_file
    self.condor_log_path = condor_log_path
    self.log_file = log_file     
    self.config_file = config_file

    self.analysis_dir = self.main_dir+'/GRB'+self.name

  # -----------------------------------------------------
  def set_addresses(self, addresses):
    self.addresses = addresses

  # -----------------------------------------------------
  def get_basic_dagname(self):
    return  self.analysis_dir+'/'+self.inifile[:-4]

  # -----------------------------------------------------
  def make_links(self, sourcedir, destdir, list_exec):
    """
    Using a list of executables names, make symbolic links rather
    than copying the files itself. That way the origin, i.e. the 
    version of the executables can be verified, at least somewhat...
    """
    
    cmd = 'cd %s;' % destdir
    for execfile in list_exec:
      cmd += 'ln -s %s/%s .;' % (sourcedir, execfile)
    system_call(self.name,cmd)


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
    self.trigger_file = 'grb%s.xml' % self.name
    utils.write_filename(xmldoc, self.trigger_file)

  # -----------------------------------------------------
  def create_call_datafind_online(self, starttime, endtime, ifo, output_location):
    """
    Creates a call to the function 'lalapps_online_datafind' to find
    the data. To be used only for data from the last 2 weeks...
    """
    executable = self.lalapps_dir+'/bin/lalapps_online_datafind'

    #type = self.type_online[ifo]
    cmd = "%s --ifo %s --gps-start-time %d --gps-end-time %d --output %s" % \
            (executable, ifo, starttime, endtime, output_location)
    return cmd


  # -----------------------------------------------------
  def create_call_datafind_offline(self,starttime, endtime, ifo, output_location):
    """
    Creates a call to the function 'ligo_data_find' to find
    the data after more than ~2 weeks. Don't ask me...
    """
    executable = self.glue_dir+'/bin/ligo_data_find --url-type file --lal-cache'

    cmd = "%s --type %s --observatory %s --gps-start-time %d --gps-end-time %d > %s" %\
          (executable, self.type_offline[ifo], ifo[0].upper(), starttime, endtime, output_location)
    return cmd

  # -----------------------------------------------------
  def run_datafind(self):
    """
    Run the datafind command to find the data
    """ 
    # create the datafind directory 
    cache_dir = "%s/GRB%s/datafind/cache" % (self.analysis_dir, self.name)
    cmd = 'mkdir -p '+cache_dir
    system_call(self.name, cmd)

    # get the start and end-time
    starttime = self.offsource_segment[0]
    endtime = self.offsource_segment[1]

    # and run the datafind command for each IFO, putting
    # the cache files directly into them
    for ifo in ifo_list:

      # create common cache-file names
      output_location = '%s/%s-DATA-%9d-%9d.cache' % (cache_dir, ifo[0].upper(), starttime, endtime)

      # decide: should I use online data (deleted after a month or so)
      # or do I require offline data? That changes everything...
      if self.use_offline_data:
        cmd = self.create_call_datafind_offline(starttime, endtime, ifo, output_location)
      else:
        cmd = self.create_call_datafind_online(starttime, endtime, ifo, output_location)

      system_call(self.name, cmd, False)

  # -----------------------------------------------------
  def check_data_to_use(self):
    """
    Checking the difference between now and the requested data
    to indentify if we can use online data or have to use
    offline data
    """ 

    # get the start time
    starttime = self.offsource_segment[0]

    # calculate the difference and test
    timediff = time.time() - offset_gps_to_linux - starttime
    self.use_offline_data = False
    if timediff>1000000:
      self.use_offline_data = True

  # -----------------------------------------------------
  def update_inifile(self, list_replacements):
    """
    Function to conveniently replace some of the parameters
    in the ini-file by specific values.
    """ 

    # read the config ini-file
    config_file = '%s/%s' % (self.analysis_dir, self.inifile)
    pc = ConfigParser.ConfigParser()
    pc.read(config_file)
     
    # Replacement loop
    for replacement in list_replacements:
      pc.set(replacement[0], replacement[1], replacement[2])

    # write out new ini-file
    cp_file = open(config_file, 'w')
    pc.write(cp_file)
    cp_file.close()


  # -----------------------------------------------------
  def prepare_analysis_directory(self):
    #
    # Now create directories and copy a bunch of files
    #

    # check if to use online or offline data
    self.check_data_to_use()
    
    # Create the main directory
    if False and os.path.exists(self.analysis_dir):
      raise IOError, "The directory %s already exists. Please (re)move"\
            " this directory or choose another name for the "\
            "analysis directory" % self.analysis_dir
    cmd = 'mkdir %s' % self.analysis_dir
    system_call(self.name, cmd)

    # copy some basic files into this directory
    cmd = 'cp %s/*.ini %s/' % \
         (self.input_dir, self.analysis_dir)
    system_call(self.name, cmd)
    cmd = 'cp %s %s/' % (self.trigger_file, self.analysis_dir)
    system_call(self.name, cmd)

    # Make some neccessary replacements in the config (ini) file
    list_replacements = [['pipeline', 'user-tag', 'GRB%s'%self.name]]
    if self.use_offline_data:
      list_replacements.append(['input', 'ligo-channel', 'LDAS-STRAIN'])
    self.update_inifile(list_replacements)

    # copy executables
    cmd = 'cd %s/bin; cp lalapps_coherent_inspiral lalapps_coherentbank \
      lalapps_coire lalapps_frjoin lalapps_inca lalapps_inspinj lalapps_inspiral \
      lalapps_inspiral_hipe lalapps_sire lalapps_thinca lalapps_tmpltbank \
      lalapps_trigbank lalapps_plot_hipe lalapps_trigger_hipe %s' %\
    (self.lalapps_dir, self.analysis_dir)
    system_call(self.name, cmd)

    #self.make_links(self, self.glue_dir+'/bin', self.analysis_dir, ['ligo_data_find','ligolw_add'])
    cmd = 'cd %s/bin; cp ligo_data_find ligolw_add %s' % \
      (self.glue_dir, self.analysis_dir)
    system_call(self.name, cmd)

    cmd = 'cd %s/bin; cp pylal_relic pylal_exttrig_llsummary plotinspiral plotnumtemplates plotthinca pylal_grbtimeslide_stats '\
          ' plotinjnum plotinspmissed plotinspinj plotsnrchi plotethinca %s' % \
      (self.pylal_dir, self.analysis_dir)
    system_call(self.name, cmd)


    # copy the segment files
    cmd = 'mv %s/*-science_grb%s.txt %s' % (self.main_dir, self.name, self.analysis_dir)
    system_call(self.name, cmd)

    #
    # make the call to trigger_hipe
    #
    cmd = 'cd %s;' % self.analysis_dir
    cmd += template_trigger_hipe % \
           (self.name, self.name, self.name, self.trigger_file, self.name, self.inifile, self.condor_log_path)
    system_call(self.name, cmd)

    #
    # similar call to set up the injection DAG
    # 
    #cmd = 'cd %s;' % self.analysis_dir
    #cmd += template_trigger_hipe_inj % \
    #       (self.name, self.name, self.name, self.trigger_file, self.name, self.inifile, self.injfile, self.condor_log_path)
    #system_call(self.name, cmd)


    # Need to unify the two cache files
    #cmd = 'cd %s/GRB%s; cat GRB%s_inj.cache >> GRB%s.cache' % \
    #  (self.analysis_dir, self.name, self.name, self.name)
    #system_call(self.name, cmd)

    # Call a subfunction to run the datafind command
    self.run_datafind()

    # Prepare the postprocesing directory at this stage and 
    # copy all files named post* into the postprocessing directory
    path = "%s/GRB%s/postprocessing" % (self.analysis_dir, self.name)
    system_call(self.name, 'mkdir -p %s/logs'%path)
    cmd = 'cp %s/post* %s' % (self.input_dir, path)
    system_call(self.name, cmd)

    # doing the same for the 'likelihood' directory 
    #path = "%s/GRB%s/likelihood" % (self.analysis_dir, self.name)
    #system_call(self.name, 'mkdir -p %s/logs'%path)
    #cmd = 'cp %s/lik* %s' % (self.input_dir, path)
    #system_call(self.name, cmd)

    # link the executables directory directly from pylal
    #exec_path = path+'/executables'
    #cmd = 'cd %s; ln -s %s/bin executables' % (path, self.pylal_dir)
    #system_call(self.name, cmd)

    # update the two DAG instances
    self.dag['onoff'] = AnalysisDag(self.name, 'onoff', self.analysis_dir)
    self.dag['inj'] = AnalysisDag(self.name, 'inj', self.analysis_dir)

    dagfile = self.get_basic_dagname()+'_onoff_uberdag.dag'
    self.dag['onoff'].set_dagname(dagfile)
    dagfile = self.get_basic_dagname()+'_inj_uberdag.dag'
    self.dag['inj'].set_dagname(dagfile)

    # create the sed file
    self.create_sed_file()

  # -----------------------------------------------------
  def prepare_injection_analysis(self):

    #
    # similar call to set up the injection DAG
    #
    cmd = 'cd %s;' % self.analysis_dir
    cmd += template_trigger_hipe_inj % \
           (self.name, self.name, self.name, self.trigger_file, self.name, self.inifile, self.injfile, self.condor_log_path)
    system_call(self.name, cmd)

    # Need to unify the two cache files
    cmd = 'cd %s/GRB%s; cat GRB%s_inj.cache >> GRB%s.cache' % \
      (self.analysis_dir, self.name, self.name, self.name)
    system_call(self.name, cmd)

    # doing the same for the 'likelihood' directory
    path = "%s/GRB%s/likelihood" % (self.analysis_dir, self.name)
    system_call(self.name, 'mkdir -p %s/logs'%path)
    cmd = 'cp %s/lik* %s' % (self.input_dir, path)
    system_call(self.name, cmd)

    # link the executables directory directly from pylal
    exec_path = path+'/executables'
    cmd = 'cd %s; ln -s %s/bin executables' % (path, self.pylal_dir)
    system_call(self.name, cmd)


  # -----------------------------------------------------
  def check_analysis_directory(self, dag_key):
    """
    Check if the dagman.out file does exist
    after some while for the dag with 'dag_key'
    """

    dag = self.dag[dag_key]

    # try it for ten times ten seconds
    for i in range(10):
      time.sleep(10)
      success = os.path.exists(dag.get_outname())
      if success:
        break

    # check the status at the end
    if not success:
      
      # send an email about this problem
      subject = 'Problems starting condor DAG'     
      email_msg = 'The condor DAG %s was not started.\n' % dag.get_dagname()
      send_mail(subject, email_msg)  

      # set the status
      if dag.status>0:
        dag.status = -dag.status
      
      return -1 
    else:

      return 1
      
  # -----------------------------------------------------
  def get_minimum_scienceseg_length(self):
    """
    Calculate the minimum science segment that
    can be used with the data given in the actual ini-file.
    The procedure below is from trigger_hipe, after
    having replaced 'cp by 'pc'
    """
    
    # get the name of the ini-file to be used
    # note: must be the inifile from CVS, just to create some information
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

      segxmlfile = "%s/segments%s_grb%s.xml" % (self.main_dir, ifo, self.name)
      segtxtfile = "%s/%s-science_grb%s.txt" % (self.main_dir, ifo, self.name)

      if not os.path.exists(segxmlfile):
        cmd = "%s/bin/ligolw_segment_query --database --query-segments --include-segments '%s' --gps-start-time %d --gps-end-time %d "\
                    "> %s" %\
                    (self.glue_dir, seg, starttime, endtime, segxmlfile)
        system_call(self.name, cmd, False)

      # 'convert' the data from the xml format to a useable format...
      # TODO: change the other places to accept the xml format
      try:
        doc = utils.load_filename(segxmlfile)
      except:
        raise IOError, "Error reading file ", segxmlfile
      segs = table.get_table(doc, "segment")
      seglist = segments.segmentlist(segments.segment(s.start_time, s.end_time) for s in segs)
      segmentsUtils.tosegwizard(file(segtxtfile, 'w'), seglist, header = True)
  
  # -----------------------------------------------------
  def get_segment_info(self,plot_segments_file = None):

    minsciseg = self.get_minimum_scienceseg_length()
    segdict = segments.segmentlistdict()

    # get the segment dicts
    for ifo in ifo_list:
      ifo_segfile = '%s/%s-science_grb%s.txt' % (self.main_dir, ifo, self.name)
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
      if offSourceSegment:
        plot.set_window(offSourceSegment, plot_offset)
      plot.highlight_segment(onSourceSegment)
      plot.finalize()
      plot.ax.set_title('Segments for GRB '+self.name)
      plot.savefig(plot_segments_file)
      plot.close()

    # store the results, cannot pickle segments objects
    if offSourceSegment:
      self.offsource_segment = [offSourceSegment[0], offSourceSegment[1]]
      self.starttime = offSourceSegment[0]
      self.endtime = offSourceSegment[1]
    self.ifolist = grb_ifolist
    if len(self.ifolist)>=2:
      self.has_data = True
    else:
      self.has_data = False

    self.ifos = ''.join(self.ifolist) 

  # --------------------------------------
  def get_html_status(self):
    """
    Returns the status of the DAGs of this instance
    in form of a dictionary.
    """

    status_dict = {}

    # loop over the dags in this instance
    for key, dag in self.dag.iteritems():
      if self.has_data:
        status = dag.get_status()
        text = dag.get_stage_name()
        if status<0:
          text = '<font color="#FF0000">%s</font>'%text
        if status==5:
          text = '<font color="#00aa00">%s</font>'%text
      else:
        text = '<font color="#666666">NoData</font>'

      # set the text
      status_dict[key]=text
       
    return status_dict

  # -----------------------------------------------------
  def calculate_optimality(self):

    # Calculate the antenna factors (for informational purpose only)
    for ifo in ifo_list:
      _, _, _, q = antenna.response(self.time, self.ra, self.de, \
                                    0.0, 0.0, 'degree', ifo)
      self.qvalues[ifo]=q

  # -----------------------------------------------------
  def get_sed_filename(self):
    """
    Returns the name of the sed file
    """
    return self.analysis_dir+'/sed.file'

  # -----------------------------------------------------
  def apply_sed_file(self, infile, outfile):
    """
    Applies the sed file to an in file
    """

    # get the sed filename
    sedfile = self.get_sed_filename()

    # run the sed command
    cmd = 'sed -f %s %s > %s' % (sedfile, infile, outfile)
    system_call(self.name, cmd, False)

  # -----------------------------------------------------
  def create_sed_file(self):
    """
    Creates the replacement sed file that will be used later
    on several in files.
    """

    # get the publishing  paths for this GRB
    publishing_path = cp.get('paths','publishing_path')
    html_path = "%s/GRB%s" % (publishing_path, self.name)


    # replace the in-file and create the DAG file
    sedfile = self.get_sed_filename()
    f = file(sedfile,'w')
    f.write("s/@GRBNAME@/GRB%s/g\n"%self.name)
    f.write("s=@ANALYSISPATH@=%s=g\n"%self.analysis_dir)
    f.write("s/@STARTTIME@/%d/g\n"%self.starttime)
    f.write("s/@ENDTIME@/%d/g\n"%self.endtime)
    f.write("s/@IFOS@/%s/g\n"%self.ifos)
    f.write("s=@LOGPATH@=%s=g\n"%self.condor_log_path)
    f.write("s/@TRIGGERTIME@/%d/g\n"%int(self.time))
    f.write("s/@RIGHTASCENSION@/%f/g\n"%float(self.ra))
    f.write("s/@DECLINATION@/%f/g\n"%float(self.de))
    f.write("s=@OUTPUTPATH@=html=g\n")
    f.write("s=@HTMLOUTPUT@=%s=g\n"%html_path)
    f.write("s/@LOGNAME@/%s/g\n" % os.getenv("LOGNAME"))
    f.write("s/@BOUNDARIESMC@/%s/g\n" % cp.get('data','mc_boundaries'))
    f.write("s/@GRBID@/%s/g\n"%self.name)
    f.write("s/@GRBPICKLE@/%s/g\n"%get_monitor_filename())
    f.write("s/@CONFIGFILE@/%s/g\n"%self.config_file)
    f.close()

