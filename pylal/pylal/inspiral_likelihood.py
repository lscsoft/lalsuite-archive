from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ilwd
from glue.ligolw import types
from glue import segments
from glue import iterutils
from pylal import grbsummary
from scipy import *
from pylal import rate
import numpy
from pylal import SearchSummaryUtils
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
    ifo_triggers[:] = [row for row in self if row.ifo == ifo]
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
    doc = utils.load_url(file, verbose = False, gz = file[-3:] == ".gz", xmldoc = ligolw.Document())
    try:
      SubTable = table.get_table(doc, InspiralLikelihoodTable.tableName)
    except:
      SubTable = []
    
    FullTable.extend(SubTable)
  return FullTable


def add_likelihood(coinc_table, likelihood_table):
  """
  Assigns likelihood to coincs using values that were stored as InspiralLikelihood table
  """
  #check if the length of coincs is equal to the length of likelihood table
  if not len(coinc_table)==len(likelihood_table):
    raise ValueError, "Error: number of coincident events does not match the number of entries in the likelihood table."\
	    "Can not combined InspiralCoinc and InspiralLikelihood tables."
    sys.exit(1)
  for (coinc, likelihood_row) in zip(coinc_table, likelihood_table):
    #check if event ids are the same
    if not coinc.event_id == likelihood_row.event_id: 
      raise ValueError, "Event ID of an element of the likelihood table does not match coinc's event ID"\
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


def loudest_coincs_in_mchirp_bins(coincs, bins):
  """Returns list  containg a coinc with maximum value of statistic in each mchirp bin.
  """ 
  loudest_coincs = []
  loudest_coincs.extend([None] * len(bins))
  num_discarded = 0
  for coinc in coincs:
    try:
      ind = bins[grbsummary.get_mean_mchirp(coinc)]
    except IndexError:
      num_discarded += 1
      continue
    if not loudest_coincs[ind] == None:
      if coinc.stat >= loudest_coincs[ind].stat:
        loudest_coincs[ind] = coinc
    else:
      loudest_coincs[ind] = coinc

  if num_discarded > 0:
        print >>sys.stderr, "warning: %d coincs fell outside "\
          "bins" % num_discarded

  return loudest_coincs


def loudest_stats_in_subsegments(coincs, full_segment, subsegment_length, ifo_combos, mchirp_bins, num_slides=None):
  """
  Return the data dictionary whose keys are ifocombos and the values are two-dimenisonal
  arrays (mchirp_bins, sebsegments) containg the loudest stat found in each subsegment and each mchirp bin.
  """
  # calculate the number of subsegments to be an integer part of division of the full segment by a subsegment.
  # This may result in the subsegments being slightly longer than the required subsegment_lenth.
  # This is done for simplicity in expectation that errors due to that are not going to be substantial.
  start_time = full_segment[0].__float__()
  stop_time = full_segment[1].__float__()  
  number_of_subsegments = int( (stop_time - start_time) // subsegment_length )
  
  #generate bins coresponding the subsegment. The last subsegment is not included, because it can be of different length.
  subseg_bins = rate.LinearBins( full_segment[0], full_segment[1], number_of_subsegments) 
  
  # Define dictionary to hold max stat arrays for each type of coincs (ifo_combo)
  max_stat_per_ifo_combo = {}
  if num_slides:
    total_number_of_trials = 2 * num_slides * number_of_subsegments
  else:
    total_number_of_trials = number_of_subsegments
  for ifo_combo in ifo_combos:
	max_stat_per_ifo_combo[combo2str(ifo_combo)] = numpy.zeros((len(mchirp_bins), total_number_of_trials))
  num_discarded_mchirp = 0
  num_discarded_subseg = 0
  for coinc in coincs:
    ifos, ifolist = coinc.get_ifos()
    try:
      mchirp_index = mchirp_bins[grbsummary.get_mean_mchirp(coinc)]
    except IndexError:
      num_discarded_mchirp += 1
    try:
      subseg_index = subseg_bins[coinc.get_time()]
    except IndexError:
      num_discarded_subseg += 1

    if num_slides:
      slide_number = coinc.slide_num
      if slide_number > 0 :
        slide_index = (slide_number - 1) * number_of_subsegments
      else:
        slide_index = ( abs(slide_number) + abs(num_slides) - 1 ) * number_of_subsegments
      max_stat_per_ifo_combo[ifos][mchirp_index, (subseg_index + slide_index)] = max(max_stat_per_ifo_combo[ifos][mchirp_index, (subseg_index + slide_index)], coinc.stat)  
    else:
      max_stat_per_ifo_combo[ifos][mchirp_index, subseg_index] = max(max_stat_per_ifo_combo[ifos][mchirp_index, subseg_index], coinc.stat)

  if num_discarded_mchirp > 0:
	  print >>sys.stderr, "warning: %d coincs fell outside "\
		" mchirp bins" % num_discarded_mchirp
  if num_discarded_subseg > 0:
	print >>sys.stderr, "warning: %d coincs fell outside "\
	  " subsegment bins" % num_discarded_subseg

	
  return max_stat_per_ifo_combo



def loudest_stat_in_slides_by_mc(slidescoincs, num_slides, mc_bins):
  """ Returns two-dimensional array (number of mc_bins, 2*num_slides) with max statistic in each 
      slide and mass chirp bin
  """
  
  loudest_stat_by_mc = numpy.zeros((len(mc_bins), 2*num_slides))
  # loop over slides
  for slide in range(1, num_slides + 1):
    forward_slide_coincs = slidescoincs.getslide(slide)
    loudest_stat_by_mc[:, (slide - 1)] = loudest_stat_in_1D_bins(forward_slide_coincs, mc_bins)	  
    backward_slide_coincs = slidescoincs.getslide(-slide)
    loudest_stat_by_mc[:, (num_slides + slide - 1)] = loudest_stat_in_1D_bins(backward_slide_coincs, mc_bins)
	  
  return loudest_stat_by_mc



def update_itemfreq(items_freq, items):
  """
  Return combined items frequency array.
  """
  new_items_freq = stats.itemfreq(items)
  if len(new_items_freq[:,0]) <= len(items_freq[:,0]):
    items_freq_1 = new_items_freq
    items_freq_2 = items_freq
  else:
    items_freq_2 = new_items_freq
    items_freq_1 = items_freq

  combined_items_freq = numpy.zeros((1,2))
  combined_items_freq[0] = items_freq_1[0]
  combined_index = 0
  index_2 = 0
  max_index_2 = len(items_freq_2[:,0]) - 1
  break_the_loop = False
  for (i,item) in enumerate(items_freq_1[:, 0]):
    while not break_the_loop:
      if index_2 <= max_index_2:
        if items_freq_2[index_2, 0] < item:
         if i > 0:
           combined_index += 1
           combined_items_freq = numpy.insert(combined_items_freq, combined_index, items_freq_2[index_2], axis=0)
         else:
           combined_items_freq = numpy.insert(combined_items_freq, combined_index, items_freq_2[index_2], axis=0)
           combined_index += 1 
         index_2 += 1
         break_the_loop = False
        elif items_freq_2[index_2, 0] == item:
          if i > 0:
            combined_index += 1
            combined_items_freq = numpy.insert(combined_items_freq, combined_index, items_freq_1[i], axis = 0)
            combined_items_freq[combined_index, 1] += items_freq_2[index_2, 1]
            index_2 += 1
          else:
            combined_items_freq[combined_index, 1] += items_freq_2[index_2, 1]
            index_2 +=1

          break_the_loop = True

        else:
          if i > 0:
            combined_index += 1
            combined_items_freq = numpy.insert(combined_items_freq, combined_index, items_freq_1[i], axis = 0)
          break_the_loop = True
      else:
        combined_index += 1
        combined_items_freq = numpy.insert(combined_items_freq, combined_index, items_freq_1[i], axis = 0)    
        break_the_loop = True
  if index_2 <= max_index_2:
     combined_items_freq = numpy.concatenate((combined_items_freq, items_freq_2[index_2:]), axis = 0)     
   
  return combined_items_freq         
         

def get_cum_prob_from_freq(freq_array):

  p1 = numpy.fliplr(freq_array.transpose())
  p1[1] = p1[1].cumsum() / numpy.sum(p1[1])
  p = numpy.fliplr(p1)
  return p



  
def get_cum_prob(x):
  """Return the array of cumulative probabilities calculated based on frequencies with which each uniq element of data array x is appearing.
     @param x: array containing values of stochastic variable X
	 The array of probabilities is a two-dimesional array. Zeroth row contains all unique values of variable X sorted in ascending order,
	 whereas first row contains the corresponding to these values cumulative probabilities. Being cumulative (and right sided or "more than")
	 the probability for the first (smallest value of X) entry must be 1. 
  """
  if len(x) >0:
    p1 = numpy.fliplr(stats.itemfreq(x).transpose())
    p1[1] = p1[1].cumsum() / numpy.sum(p1[1])
    p = numpy.fliplr(p1)
    return p
  else:
     raise ValueError, "empty frequency array"

  
  
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
    if len(stat_list) > 0:
      prob_list.append(get_cum_prob(numpy.asarray(stat_list)))
    else: 
      prob_list.append(None)
  return prob_list	
	  
  	 
# Convert ifo list to string
def combo2str(combo):
  ifo_list = ""
  for ifo in combo:
    ifo_list += ifo
  return ifo_list


def get_segments_analyzed_from_file(file):
  """
  Return list of segments that were analyzed to generate the file.
  @param file: name (including location) of ligolw xml file 
  """
  segDict = SearchSummaryUtils.GetSegListFromSearchSummaries([file])
  seglist = segments.segmentlist(iterutils.flatten(segDict.values()))

  return seglist


def weighted_histogram(array,weight=None, norm=False, nbins=10, set_max=None, set_min=None):
  """
  Weighted histogram function; It returns 
  y is a numpy array containing number of events in the corresponding bin;
  if norm=True then y is normalized to unity.
  x is a numpy array containing bins (namely center of each bin).
  @param array: numpy array containing the quatities for which the histogram will be generated.
  @param bins: integer that sets the number of bins in the histogram.
  @param weight: if given must be array of weights for each element of array to be binned.
  @param norm: if True then the histogram is normed to 1, otherwise nopthing is done.
  @param set_max: parameter that sets the maximum value for bins, if not given it is determined from the data.
  @param set_min: parameter that sets the minimum value for bins, if not given it is determined from the data.
  """
  
  if weight==None:
    weight = numpy.ones(len(array))
  if not set_max:	
    max_element = numpy.max(array)
  else:
    max_element = set_max
  if not set_min:
    min_element = numpy.min(array)
  else:
    min_element = set_min
  bins_object = rate.LinearBins(min_element, max_element, nbins)
  y = numpy.zeros(nbins)
  for i in range(len(array)):
    if (array[i] >= min_element) and (array[i] <= max_element):
      bin_index = bins_object[array[i]]
      y[bin_index] += weight[i]
  if norm:
    normalization = numpy.sum(y)
    y = y/normalization
  x = bins_object.centres()
  return y, x

