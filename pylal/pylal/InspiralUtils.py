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
import sys
import copy

from glue.ligolw import utils
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import SnglInspiralUtils

# set default color code for inspiral plotting functions
colors = {'G1':'k','H1':'r','H2':'b','L1':'g','V1':'m'}
symbols = {'G1':'Y','H1':'x','H2':'o','L1':'+','V1':'1'}

class InspiralPage(object):
  """
  This is a class to contain all the bits of a inspiral page
  showing the results of a piece of code.
  """

  def __init__(self, options):
    """
    Initializes this class with the options.
    """
    self.opts = options
    
    self.fname_list = []
    self.tag_list = []
    self.html_footer = ""

    # just adding some stuff to the opts structure
    # (should be fixed later)
    initialise(self.opts, os.path.basename(sys.argv[0]))

  def add_plot(self, plot_fig, text):
    """
    Add a plot to the page
    """
    
    fname = set_figure_name(self.opts, text)
    fname_thumb = savefig_pylal(fname, fig=plot_fig)
    
    self.fname_list.append(fname)
    self.tag_list.append(fname)

  def write_page(self):
    """
    create the page
    """
    if self.opts.enable_output:
      html_filename = write_html_output(self.opts, sys.argv[1:],\
                                        self.fname_list, self.tag_list,\
                                        comment=self.html_footer or None)
      write_cache_output(self.opts, html_filename, self.fname_list)

  def write(self, text):
    """
    Write some text to the standard output AND
    to the page.
    """
    print text
    self.html_footer+=text+'<br>'
      
  

def savefig_pylal(filename=None, filename_thumb=None, doThumb=True, dpi=None,
  dpi_thumb=50, fig=None):
  """
  @param filename: filename in which to save the figure
  @param filename_thumb: filename into which save a thumbnail of the figure
  @param doThumb: save the thumbnail or not (doThumb=True by default)
  @param dpi: resolution of the figure
  @param dpi_thumb: resolution of the thumbnail (dpi=50 by default)
  @param fig: the particular figure you wish to save (current figure by
              default)
  @return filename_thumb if a thumbnail was created (computed from filename
          by default)

  """
  import pylab
  
  # fill in non-trivial defaults
  if fig is None:
    fig = pylab.gcf()
  if dpi is None:
    dpi = pylab.rcParams["savefig.dpi"]
  if doThumb and (filename_thumb is None):
    if filename is None:
      raise ValueError, "must provide filename_thumb or filename if doThumb "\
        "is True"
    index = filename.rindex('.')
    filename_thumb = filename[0:index] + '_thumb' + filename[index:]
  
  # save picture into a file
  if filename is not None:
    fig.savefig(filename, dpi=dpi)

  # save thumbnail into a file if requested
  if doThumb:
    fig.savefig(filename_thumb, dpi=dpi_thumb)

  return filename_thumb


def ErrorMessagePlotting(opts, thisplot):
   """

   """
   text = "---Error in "+opts.name+"in plotting functions "+thisplot
   if "chi" in thisplot:
     text += "\n---possible reasons related to chi-square (are you reading first stage triggers ?)"
   print >> sys.stderr, text


def message(opts, text):
  """

  """
  if opts.verbose:
    print text
  return text+'</br>\n'

def set_figure_name(opts, text):
  """
  return a string containing a standard output name for pylal 
  plotting functions.
  """
  fname = "Images/" + opts.prefix + "_"+text + opts.suffix + ".png"
  
  if opts.output_path is not None:
    fname = opts.output_path + fname

  return fname

def write_html_output(opts, args, fnameList, tagLists, \
			doThumb=True, cbcweb = False, mapList = [],\
			comment=None ):
  """
  @param opts: The options from the calling code
  @param args: The args from the calling code
  @param fnameList: A list of the filenames
  @param tagLists: A list for the tags, getting added to the links
  @param doThumb: Uses the _thumb file as the sourcs for the images
  @param cbcweb: Creates the output as a CBC webpage
  @param mapList: A list of dictionaries to create the image maps
  """

  # -- the HTML document and output cache file
  # -- initialise the web page calling init_page
  page, extra = init_markup_page(opts)
  if cbcweb:
    page.addheader("<%method title>" + opts.name + " results</%method>")
    page.addheader("<%method headline>" + opts.name + " results</%method>")
    page.addheader("<%method cvsid> $Id$ </%method>")
  else:
    page.h1(opts.name + " results")

  page.p(opts.prefix + opts.suffix)
  page.hr()

  # -- filename
  if cbcweb:
    html_filename = opts.prefix + opts.suffix +"_publish.html"
  else:
    html_filename = opts.prefix + opts.suffix +".html"  
  if opts.output_path:
    html_filename = opts.output_path + html_filename
  html_file = file(html_filename, "w")

  # loop over the contents
  for tag,filename in zip(tagLists,fnameList):

    # set the correct name for linking (two '//' does not bother)
    if cbcweb:
      fname = opts.html_for_cbcweb + "/Images/" + os.path.basename(filename)
    else:
      fname = "Images/" + os.path.basename(filename)
     

      # set the thumbnail pictures if required
    if doThumb:
      fname_thumb = fname[:-4] + "_thumb.png"
    else:
      fname_thumb =fname

    # add the image to tge page
    page.a(extra.img(src=[fname_thumb], width=400, \
        alt=tag, border="2"), title=tag, href=[ fname])
    
  page.add("<hr/>")

  # add maps to this page
  if len(mapList)>0:
    m=0
    for mapDict in mapList:
      m+=1
      page.add( mapDict['text']+'<br>' )
      page.add( '<IMG src="%s" width=800px '\
                'usemap="#map%d">' % ( mapDict['object'], m) )
      page.add( '<MAP name="map%d"> <P>' % m )
      n=0
      for px, py, link in zip( mapDict['xCoords'],  \
                               mapDict['yCoords'],  \
                               mapDict['links']):
        n+=1
        page.add( '<area href="%s" shape="circle" '\
                  'coords="%d, %d, 5"> Point%d</a>' %\
                  ( link, px, py, n) )
      page.add('</P></MAP></OBJECT><br>')
      page.add("<hr/>")    

  if opts.enable_output:
    if comment is not None:
      page.add("<div> "+comment+"</div>")
      page.hr()
    if CoincSummTable is not None:
      page.add(CoincSummTable)
      page.hr()
    text = writeProcessParams( opts.name, opts.version,  args)
    page.add(text)
    if comment is not None:
      page.hr()
      page.add("<div> "+comment+"</div>")
    html_file.write(page(False))
    html_file.close()

  return html_filename


def write_cache_output(opts, html_filename,fnameList):
  """
  write the output cache file of theplotting functions
  """

  output_cache_name = opts.prefix + opts.suffix +'.cache'
  if opts.output_path:
    output_cache_name = opts.output_path + output_cache_name
  this = open(output_cache_name, 'w')
  if opts.enable_output:
    this.write(os.path.basename(html_filename) + '\n')
  for filename in fnameList:
    if str(filename).endswith('.png'): 
      fname = "Images/"+os.path.basename(filename) # set the correct name for linking
    elif str(filename).endswith('.html'): 
      fname = os.path.basename(str(filename)) # set the correct name for linking
    this.write(fname + '\n')
  this.close()


def writeProcessParams(name, version, command): 
  """
  Convert input parameters from the process params that the code was called 
  with into a formatted string that can be saved within an other document 
  (e.g., HTML)

  @param name: name of the executable/script
  @param version:version of the executable/script
  @param command: command line arguments from a pylal script
  @return text
  """
  text = "Figure(s) produced with " + name + ", " \
      + version + ", invoked with the following command line arguments:" \
      + '</br>\n<p style="width:80%; color:blue">'+ name
  for arg in command:
    text += " " +  arg
  text+='</p>'
  
  return text

def AddFileToCache(fname, cache):
  """
  Add the given file to the lal.Cache

  @param fname:
  @param cache:
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

  @param fileList : a list of file
  @return cache
  """
  cache = lal.Cache()
  for file in fileList:
    AddFileToCache(file, cache)
  return(cache)


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



def initialise(opts, name, version):
  """
  Create suffix and prefix that will be used to name the output files.

  @param opts : the user arguments (user_tag, gps_end_time and 
  gps_start_time are used).
  @param name: name of the calling function/executable
  @return prefix 
  @return suffix
  """

  # compose prefix
  prefix = name
  try:
    if opts.ifo_times:
      prefix = opts.ifo_times +"-"+ prefix
  except:
     print >> sys.stderr, "--ifo-time option not implemented in the "+name +" executable. skipping..."
     pass
  try:
    if opts.ifo_tag:
      prefix = prefix + "_" + opts.ifo_tag
  except: 
     print >> sys.stderr, "--ifo-tag option not implemented in the "+name +" executable. skipping..."
     pass
  try:
    if opts.user_tag:
      prefix = prefix + "_" + opts.user_tag
  except: 
     print >> sys.stderr, "--user-tag option not implemented in the "+name +" executable. skipping..."
     pass
  

  # compose suffix
  try:
    if opts.gps_start_time and opts.gps_end_time :
      suffix = "-"+str(opts.gps_start_time)+"-"+str(opts.gps_end_time-opts.gps_start_time)
    else:
      suffix = "-unspecified-gpstime"
  except:
     suffix = "-unspecified-gpstime"
     print >> sys.stderr, "--gps-start-time and/or --gps-end-time  option not implemented in the "+\
           name +" executable. skipping..."
     pass
  

  opts.prefix = prefix
  opts.suffix = suffix
  opts.name = name
  opts.version = version

  # make sure output_path is set correctly
  if opts.output_path is not None:
    opts.output_path = opts.output_path +'/'

    # create output file if required
    if not os.path.exists( opts.output_path ):
      os.mkdir (opts.output_path)

    if not os.path.exists( opts.output_path+"Images" ):
      os.mkdir (opts.output_path+"Images")
      
  else:
    if not os.path.exists( "Images" ):
      os.mkdir( "Images")


  
  return opts

def init_markup_page( opts):
  """
  Load the markup module, and initialise the HTML document if the opts 
  argument contains enable_ouput option.

  @param  opts : the user arguments 
  @return page 
  @return extra 
  """
  # Initialise the html output file
  if opts.enable_output:
    try:
      from glue import markup
      from glue.markup import oneliner as extra_oneliner
    except:
      raise ImportError("Require markup.py to generate the html page")

    page = markup.page()
    try:
      page.init(title=__title__)
    except:
      page.init()

  return page, extra_oneliner


def readHorizonDistanceFromSummValueTable(fList, verbose=False):
  """
  read in the SummValueTables from a list of files and return the
  horizon distance versus total mass

  @param fList:       list of input files
  @param verbose: True of False (default is False)
  """

  output = {}
  massOutput = {}
  count = 0
  if len(fList) == 0:
    return output

  # for each file in the list 
  for thisFile in fList:
    if verbose:
      print str(count+1)+"/"+str(len(fList))+" " + thisFile
    count = count+1
    massNum = 0

    doc = utils.load_filename(thisFile, gz = thisFile.endswith(".gz"))
    try:
      summ_value_table = table.get_table(doc, lsctables.SummValueTable.tableName)
    except ValueError:
      print "ValueError in readHorizonDistanceFromSummValueTable whiile reading summvalue table from file ", thisFile
      return output,massOutput

    # if not summ_value table was filled , then simply returns 
    if summ_value_table is None:
      return output,massOutput
    
    # else
    for row in summ_value_table:
      # we should fnd a name "inspiral_Effective_distance
      if row.name == 'inspiral_effective_distance':
        # it may be that the file read is an inspiral file containing only the BNS infomration
        if (row.comment == '1.40_1.40_8.00') or (row.comment == '1.4_1.4_8'):
          if not output.has_key(row.ifo):
            output[row.ifo] = lsctables.New(lsctables.SummValueTable)
          output[row.ifo].append(row)
        # or a template bank containing a whole list of inspiral_effective_distance
        else:
          if not massOutput.has_key(row.ifo):
            massOutput[row.ifo] = [lsctables.New(lsctables.SummValueTable)]
          if len(massOutput[row.ifo]) < massNum + 1:
            massOutput[row.ifo].append(lsctables.New(lsctables.SummValueTable))
          massOutput[row.ifo][massNum].append(row)
          massNum += 1
  return output,massOutput



def isPlayground(table):
  """
  @param gpsTime: a list of valid GPS time
  @return True if it belongs to the playground time
  """
  start= 729273613 ;
  interval = 6370;
  len = 600;
  gpsTime = table.geocent_end_time + table.geocent_end_time_ns * 1e-9
  if ((gpsTime-start) % interval)<len:
    return True
  else:
    return False

