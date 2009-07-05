#!/usr/bin/python

import sys
import os
import urlparse
from optparse import *
from subprocess import Popen,PIPE


from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils

##############################################################################
#
#          definition of the LVAlertTable
#
##############################################################################

class LVAlertTable(table.Table):
  """
  for reference, file is written as
  file://host/path_to_file/file
  and uid is the unique id assigned by gracedb
  """
  tableName = "LVAlert:table"
  validcolumns = {
    "file": "lstring",
    "uid": "lstring"
    }
    
class LVAlertRow(object):
  __slots__ = LVAlertTable.validcolumns.keys()
  
LVAlertTable.RowType = LVAlertRow

##############################################################################
#
#          useful utilities
#
##############################################################################

def _parse_file_url(file_url):
  """
  simple function to parse the file urls of the form:
  file://host_name/path_to_file/file
  where path_to_file is assumed to be the /private subdir of a gracedb entry
  returns:
  host: host_name in the above example
  full_path: path_to_file/file in the above example
  general_dir: the /general subdir of the gracedb entry (where data not
  produced by the event supplier should be placed)
  """
  parsed = urlparse.urlparse(file_url)
  host = parsed[1]
  full_path = parsed[2]
  general_dir = os.path.join(os.path.split(os.path.split(full_path)[0])[0], \
                             'general')
  return host, full_path, general_dir

def get_LVAdata_from_stdin(std_in):
  """
  this function takes an LVAlertTable from sys.stdin and it returns:
  host: the machine the payload file was created on
  full_path: the full path to (and including) the payload file
  general_dir: the directory in gracedb that the output of your code should
               be written to
  uid: the gracedb unique id associated with the event in the LVAlertTable
  """
  doc = utils.load_fileobj(std_in)[0]
  lvatable = table.get_table(doc, LVAlertTable.tableName)
  host, full_path, general_dir = _parse_file_url(lvatable[0].file)
  uid = lvatable[0].uid

  return host, full_path, general_dir, uid

def get_LVAdata_from_file(file):
  """
  this function takes the name of an xml file containing a single LVAlertTable
  and it returns:
  host: the machine the payload file was created on
  full_path: the full path to (and including) the payload file
  general_dir: the directory in gracedb that the output of your code should
               be written to
  uid: the gracedb unique id associated with the event in the LVAlertTable
  """
  doc = utils.load_filename(opts.online_input)
  lvatable = table.get_table(doc, LVAlertTable.tableName)
  lvatable = table.get_table(doc, LVAlertTable.tableName)
  host, full_path, general_dir = _parse_file_url(lvatable[0].file)
  uid = lvatable[0].uid

  return machine, full_path, general_dir, uid


#the following is meant as a template for small jobs
#notes:
#   *we only use the vanilla universe which is appropriate for python
#    jobs and things not condor-compiled
#   *it only works for single-process jobs anything more complicated will
#    require a dag
condor_sub_template = \
                    """
                    universe = vanilla
                    executable = macroexecutible
                    arguments = macroargs
                    log = macrolog
                    error = macroerr
                    output = macroout
                    getenv = True
                    notification = never
                    queue
                    """

def write_condor_sub(executible, args, logdir, uid):
  """
  write a simple condor submission file
  uid: unique id used in naming the files (to avoid conflicts)
  executible: the name of the executible file
  args: a list of arguments to executible
  logdir: directory to keep log files
  returns the name of the file
  """
  subfile = condor_sub_template
  subfile.replace('macroexecutible', executible)
  subfile.replace('macroargs', args)
  subfile.replace('macrolog', os.join(logdir,str(uid)+'.log'))
  subfile.replace('macroerr', os.join(logdir,str(uid)+'.err'))
  subfile.replace('macroout', os.join(logdir,str(uid)+'.out'))

  fname = str(uid) + '.sub'
  f = open(fname,'w')
  f.write(subfile)
  f.close()

  return fname

def submit_condor_job(subfile):
  """
  submit the subfile to condor
  returns the process id
  """
  p = Popen(["condor_submit "+subfile], shell=True).pid
  
  return p

