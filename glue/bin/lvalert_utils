#!/usr/bin/python

import sys
import os
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
  machine:full/path/to/file.ext
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

def get_LVAdata_from_stdin(std_in):
  """
  this function takes an LVAlertTable from sys.stdin and it returns:
  machine: the machine the payload file was created on
  full_path: the full path to (and including) the payload file
  general_dir: the directory in gracedb that the output of your code should
               be written to
  uid: the gracedb unique id associated with the event in the LVAlertTable
  """
  doc = utils.load_fileobj(std_in)[0]
  lvatable = table.get_table(doc, LVAlertTable.tableName)
  machine, full_path = lvatable[0].file.split(':')
  general_dir = os.path.join(os.path.split(os.path.split(full_path)[0])[0],'general')
  uid = lvatable[0].uid

  return machine, full_path, general_dir, uid

def get_LVAdata_from_file(file):
  """
  this function takes the name of an xml file containing a single LVAlertTable
  and it returns:
  machine: the machine the payload file was created on
  full_path: the full path to (and including) the payload file
  general_dir: the directory in gracedb that the output of your code should
               be written to
  uid: the gracedb unique id associated with the event in the LVAlertTable
  """
  doc = utils.load_filename(opts.online_input)
  lvatable = table.get_table(doc, LVAlertTable.tableName)
  machine, full_path = lvatable[0].file.split(':')
  general_dir = os.path.join(os.path.split(os.path.split(full_path)[0])[0],'general')
  uid = lvatable[0].uid

  return machine, full_path, general_dir, uid


#this stuff works for single-process jobs
#anything more complicated will require a dag
condor_sub_template = \
                    """
                    universe = standard
                    executable = macroexecutible
                    arguments = macroargs
                    log = macrolog
                    error = macroerr
                    output = macroout
                    getenv = True
                    notification = never
                    queue
                    """

def write_condor_sub(uid, executible, args, logdir):
  """
  write a simple condor submission file
  uid: unique id (probably from gracedb)
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

