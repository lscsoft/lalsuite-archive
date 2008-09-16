from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ilwd
from glue.ligolw import types
from pylal import grbsummary
from scipy import *
from pylal import rate
################################################################################
# Definitions of some of the classes and functions that are used in the main
# program.
################################################################################
class InspiralLikelihoodTable(table.Table):
  tableName = "Inspiral_Likelihood:table"
  validcolumns = {
    "ifo": "lstring",
    "false_alarm_prob": "real_4",
    "detect_prob": "real_4",
    "likelihood": "real_4",
    "bg_frac": "real_4",
    "detect_bg_frac_ratio": "real_4",
    "event_id": "ilwd:char"
  }
  def get_column(self,column):
    return self.getColumnByName(column).asarray()

  def getslide(self,slide_num):
    """
    Return the triggers with a specific slide number.
    @param slide_num: the slide number to recover (contained in the event_id)
    """
    slideTrigs = lsctables.New(InspiralLikelihoodTable)
    slideTrigs.extend([r for r in self if r.get_slide_number() == slide_num])
    return slideTrigs

  def coinctype(self, ifo):
    """
    Return the events with the given IFO.
    @param ifo: parameter specifying the IFO ( e.g. "H1H2L1")  
    """
    ifo_triggers=lsctables.New(InspiralLikelihoodTable)

    for row in self:
      if row.ifo == ifo:
        ifo_triggers.append(row)

    return ifo_triggers

class InspiralLikelihood(object):
  __slots__ = InspiralLikelihoodTable.validcolumns.keys()

  def get_id_parts(self):
    """
    Return the three pieces of the int_8s-style sngl_inspiral
    event_id.
    """
    int_event_id = int(self.event_id)
    x = int_event_id // 1000000000
    slidenum = (int_event_id % 1000000000) // 100000
    y = int_event_id % 100000
    return x, slidenum, y

  def get_slide_number(self):
    """
    Return the slide-number for this trigger
    """
    x, slidenum, y = self.get_id_parts()
    slide_number = slidenum
    if slide_number > 5000:
      slide_number = 5000 - slide_number
    return slide_number

InspiralLikelihoodTable.RowType = InspiralLikelihood


def ReadInspiralLikelihoodFromFiles(fileList):
  """
  Read the InspiralLikelihoodTables from a list of files
  @param fileList: list of input files
  """
  FullTable = lsctables.New(InspiralLikelihoodTable)
  for file in fileList:
    doc = utils.load_url(file, verbose = False, gz = file[-3:] == ".gz",
        xmldoc = ligolw.Document())
    try:
      SubTable = ligolw.table.get_table(doc, InspiralLikelihoodTable.tableName)
    except:
      SubTable = None
    if SubTable:
      for row in SubTable:
        FullTable.append(row)
  return FullTable


def add_likelihood(coinc_table, likelihood_table):
  """
  Assigns likelihood to coincs using values that were stored as InspiralLikelihood table
  """
  #check if the length of coincs is equal to the length of likelihood table
  if not len(coinc_table)==len(likelihood_table):
    print >>sys.stderr, "Error: number of coincident events does not match the number of entries in the likelihood table."\
	    "Can not combined InspiralCoinc and InspiralLikelihood tables."
    sys.exit(1)
  for (coinc, likelihood_row) in zip(coinc_table, likelihood_table):
    #check if event ids are the same
    if not coinc.event_id == likelihood_row.event_id: 
      print >>sys.stderr, "Event ID of an element of the likelihood table does not match coinc's event ID"\
	  "Can not merge InspiralCoinc and InspiralLikelihood tables"
      sys.exit(1)
    coinc.likelihood = likelihood_row.likelihood  




def generate_prefix_and_suffix(ifo_times="", program_name="", ifo_tag="",
  user_tag="", gps_start_time=0, gps_end_time=0):
  """
  Generates the string according to inspiral pipeline conventions, that can
  be used as a name of a file.
  """
  prefix = program_name
  if ifo_times:
    prefix = ifo_times + "-" + prefix
  if ifo_tag:
    prefix = prefix + "_" + ifo_tag
  if user_tag:
    prefix = prefix + "_" + user_tag

  if gps_start_time > 0 and gps_end_time > 0:
    suffix = "_" +str(gps_start_time) + "-" + \
      str(gps_end_time - gps_start_time)
  else:
    suffix = ""
  return prefix, suffix


def loudest_stat_in_1D_bins(coincs, bins):
  """Returns one-dimensional array containg maximum value of statistic in each bin.
  """ 
  loudest_stat = numpy.array(len(bins))
  num_discarded = 0
  for coinc in coincs:
    try:
      ind = mc_bins[grbsummary.get_mean_mchirp(coinc)]
    except IndexError:
      num_discarded += 1
      continue
    loudest_stat[ind] = max(coinc.stat, loudest_stat[ind])

  if num_discarded > 0:
        print >>sys.stderr, "warning: %d coincs fell outside "\
          "bins" % num_discarded

  return loudest_stat





def loudest_stat_in_slides_by_mc(slidescoincs, num_slides, mc_bins):
  """ Returns two-dimensional array (number of mc_bins, 2*num_slides) with max statistic in each 
      slide and mass chirp bin
  """
  
  loudest_stat_by_mc = numpy.zeros((len(mc_bins), 2*num_slides))
  # loop over slides
  for slide in range(1, opts.num_slides + 1):
    forward_slide_coincs = slidescoincs.getslide(slide)
    loudest_stat_by_mc[:, (slide - 1)] = loudest_stat_in_1D_bins(forward_slide_coincs, mc_bins)	  
    backward_slide_coincs = slidescoincs.getslide(-slide)
    loudest_stat_by_mc[:, (num_slides + slide - 1)] = loudest_stat_in_1D_bins(backward_slide_coincs, mc_bins)
	  
  return loudest_stat_by_mc



  
def get_cum_prob(x):
  """Return the array of cumulative probabilities calculated based on frequencies with which each uniq element of data array x is appearing.
     @param x: array containing values of stochastic variable X
	 The array of probabilities is a two-dimesional. Zeroth row contains all unique values of variable X sorted in ascending order,
	 whereas first row contains the corresponding to these values cumulative probabilities. Being cumulative (and right sided or "more than")
	 the probability for the first (smallest value of X) entry must be 1. 
  """
  p1 = numpy.fliplr(stats.itemfreq(x).transpose())
  p1[1] = p1[1].cumsum() / numpy.sum(p1[1])
  p = numpy.fliplr(p1)
  return p 
  
  
def get_prob_by_mc(loudest_stat_by_mc):
  """Return a list wich elements are two-dimensional arrays of cumulative probabilites for each of chirp mass bin.
     Length of the list is equal to number of bins.
	 @param loudest_stat_by_mc: two-dimensional array containg the loudest stat in each bin, trial
  """
  prob_list = []
  
  for i in range(len(loudest_stat_by_mc[:,0])):
    prob_list.append(get_cum_prob(loudest_stat_by_mc[i]))
  return prob_list
  
def get_prob_by_mc_for_injections(stat_lists):
  """Return a list wich elements are two-dimensional arrays of cumulative probabilites for each of chirp mass bin.
     Length of the list is equal to number of bins.
	 @param stat_lists: list of lists of injection stats in each bin.
  """		
  prob_list = []
  
  for stat_list in stat_lists:
    prob_list.append(get_cum_prob(numpy.asarray(stat_list)))
  return prob_list	
	  
  	 
# Convert ifo list to string
def combo2str(combo):
  ifo_list = ""
  for ifo in combo:
    ifo_list += ifo
  return ifo_list



		
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
