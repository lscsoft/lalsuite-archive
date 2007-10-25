#!/usr/bin/python
"""
Utilities for the inspiral plotting functions
"""
__version__ = "$Revision$"
__date__ = "$Date$"
__Id__ = "$Id$"

# $Source$

from glue import lal
from glue import segments
import socket, os

from glue.ligolw import utils

from glue.ligolw import lsctables



def writeProcessParams(name, version, command_line): 
  """
  function to write out the process params that the code was called with
  """
  text = "Figures produced with "+name+", " \
      + version[1:len(version)-1] + ", invoked with arguments:\n\n" \
      + name
  for arg in command_line:
    text += " " +  arg
  
  return text

def AddFileToCache(fname, cache):
  """
  Add the given file to the lal.Cache
  """
  file_name = fname.split('.')[0].split('-')
  cache.append(lal.CacheEntry( file_name[0], file_name[1],
    segments.segment(int(file_name[2]), 
      int(file_name[2]) + int(file_name[3])),
    'file://' + socket.gethostbyaddr(socket.gethostname())[0] + \
     os.getcwd() + '/' + fname))

def GenerateCache(fileList):
  """
  Generate a lal.Cache for the list of files
  """
  cache = lal.Cache()
  for file in fileList:
    AddFileToCache(file, cache)
  return(cache)


#def ContentHandler(ligolw.PartialLIGOLWContentHandler):
def ContentHandler(PartialLIGOLWContentHandler):
  """
  
  """
  def __init__(self, xmldoc):
    """
    New content handler that only reads in the SummValue table
    """
    def element_filter(name, attrs):
      """
      Return True if name and attrs describe a SummValueTable
      """
      return lsctables.IsTableProperties(lsctables.SummValueTable, name, attrs) 
    PartialLIGOLWContentHandler.__init__(self, xmldoc, element_filter)


def create_output_name(opts, name):
  """
  Create the suffix and prefix used for the naming convention.
  @param  opts : the user arguments (user_tag, gps_end_time 
        and gps_start_time are used)
  """
  if not opts.user_tag:
    prefix = opts.ifo_times +"-"+ name + "_"
  else:
    prefix = opts.ifo_times +"-"+ name + "_" + opts.user_tag + "_"
  if opts.gps_start_time and opts.gps_end_time :
    suffix = "-"+str(opts.gps_start_time)+"-"+str(opts.gps_end_time-opts.gps_start_time)
  else:
    suffix = "-unspecified-gpstime"

  return prefix, suffix

def init_markup_page( opts):
  """
  Load the markup module, and initialise the HTML document if the opts 
  argument contains enable_ouput option.
  @param  opts : the user arguments 
  @return page : the HTML document
  @return extra : the onliner (see markup.py documentation)
  """
  # Initialise the html output file
  if opts.enable_output is True:
    try:
      import markup
      from markup import oneliner as extra_oneliner
    except:
      raise ImportError("Require markup.py to generate the html page")

    page = markup.page()
    try:
      page.init(title=__title__)
    except:
      page.init()

  return page, extra_oneliner


def readFiles(fList, verbose=False):
  """
  read in the SimInspiralTables from a list of files
  @param fList:       list of input files
  """
  output = {}
  massOutput = {}
  count = 0
  if len(fList) == 0:
    return output

  for thisFile in fList:
    if verbose is True:
      print str(count)+"/"+str(len(fList))+" " + thisFile
    count = count+1
    massNum = 0
    doc = utils.load_filename(thisFile, gz = thisFile.endswith(".gz"))
    for row in doc.childNodes[0]:
      if row.name == 'inspiral_effective_distance':
        if (row.comment == '1.40_1.40_8.00') or (row.comment == '1.4_1.4_8'):
          if not output.has_key(row.ifo):
            output[row.ifo] = lsctables.New(lsctables.SummValueTable)
          output[row.ifo].append(row)
    for row in doc.childNodes[0]:
      if row.name == 'inspiral_effective_distance':
        if not massOutput.has_key(row.ifo):
          massOutput[row.ifo] = [lsctables.New(lsctables.SummValueTable)]
        if len(massOutput[row.ifo]) < massNum + 1:
          massOutput[row.ifo].append(lsctables.New(lsctables.SummValueTable))
        massOutput[row.ifo][massNum].append(row)
        massNum += 1



  return output,massOutput






