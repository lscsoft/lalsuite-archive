import sys
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal.inject import light_travel_time
from pylal.tools import XLALCalculateEThincaParameter
from pylal.xlal import date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
import glue.iterutils
import numpy
import cmath

########################################
# the list of available IFOs comes up several times
ifos = ("G1", "H1", "H2", "L1", "T1", "V1")

########################################
# helper functions

def get_ifo_combos(ifo_list):
  ifo_combos = []
  for num_ifos in range(2, len(ifo_list) + 1):
    ifo_combos.extend(list(glue.iterutils.choices(ifo_list, num_ifos)))

  return ifo_combos

def simpleEThinca(trigger1, trigger2):
  """ 
  Return the e-thinca corresponding to the distance in  parameter space between two inspiral triggers.

  The average distance defined below  is only an approximation to the true distance and is
  valid whenever two triggers are nearby. The simplified version of
  the e-thinca parameter is calculated based on that definition of the distance.

  d_average=(1/2)[(Gamma(x1)_{ij}(x2-x1)^i(x2-x1)^j)^(1/2) + (Gamma(x2)_{ij}(x2-x1)^i(x2-x1)^j)^(1/2)]
  then simple_ethinca= d_average^2/4  
  
  @param trigger1: is a single inspiral triggers.
  @param trigger2: is a single inspiral triggers.
  """ 
  #dend_time = (trigger2.end_time - trigger1.end_time) + (trigger2.end_time_ns - trigger1.end_time_ns)*10**(-9)

  # HACK
  #dend_time = numpy.sign(dend_time)*max(0.0, abs(dend_time) - light_travel_time(trigger1.ifo, trigger2.ifo))
  #dt = (trigger2.end_time - trigger1.end_time) + (trigger2.end_time_ns - trigger1.end_time_ns)*10**(-9)
  #FIX ME end_time for time slides is poorly defined, we should sort it out
  #dend_time = (trigger2.end_time_ns - trigger1.end_time_ns)*10**(-9)

  dt = 0
  dtau0 = trigger2.tau0-trigger1.tau0
  dtau3 = trigger2.tau3-trigger1.tau3

  dist1 = dt * (dt * trigger1.Gamma0 + dtau0 * trigger1.Gamma1 + dtau3 * trigger1.Gamma2) + \
       dtau0 * (dt * trigger1.Gamma1 + dtau0 * trigger1.Gamma3 + dtau3 * trigger1.Gamma4) + \
       dtau3 * (dt * trigger1.Gamma2 + dtau0 * trigger1.Gamma4 + dtau3 * trigger1.Gamma5)
  
  dist2 = dt * (dt * trigger2.Gamma0 + dtau0 * trigger2.Gamma1 + dtau3 * trigger2.Gamma2) + \
       dtau0 * (dt * trigger2.Gamma1 + dtau0 * trigger2.Gamma3 + dtau3 * trigger2.Gamma4) + \
       dtau3 * (dt * trigger2.Gamma2 + dtau0 * trigger2.Gamma4 + dtau3 * trigger2.Gamma5)

  average_distance = 0.5 * (numpy.sqrt(dist1) + numpy.sqrt(dist2))

  simple_ethinca = (average_distance * average_distance) / 4.0
  
  return simple_ethinca
  
def convert_to_polar_coord(x, y):  
  r_rms = ((x**2 + y**2)/2.0)**(0.5)
  phi = numpy.arctan2(float(y),float(x))
  return r_rms, phi

def convert_to_spherical_coord(x, y, z):
   """ 
   This includes rotation by 45 deg around z follwed by rotation by 45 deg
   around y to define new x,y,z coordinates in which z would be along (1,1,1) 
   vector of the original coordinates. The spherical coordinates then defined with respect
   this new set of Cartesian coordinates.
   """
   sq_root_three = 3**(0.5)
   sq_root_two = 2**(0.5)
   x_new = (x + y - 2.0*z) / (sq_root_three*sq_root_two)
   y_new = (y - x) / sq_root_two
   z_new = (x + y + z) / sq_root_three

   r_rms = ((x_new**2 + y_new**2 + z_new**2)/3.0)**(0.5)
   if (x_new != 0.0 or y_new != 0.0) and (z_new == 0.0):
     theta = cmath.pi/2.0
   else:
     theta = numpy.arctan((x_new**2 + y_new**2)**(0.5)/float(z_new))
   phi = numpy.arctan2(float(y_new),float(x_new))
   return r_rms, theta, phi

def frac_error_sq(candidate, trigger):
  return ((float(candidate) - float(trigger)) / float(candidate))**2


def readCoincInspiralFromFiles(fileList,statistic=None):
  """
  read in the Sngl and SimInspiralTables from a list of files
  if Sngls are found, construct coincs, add injections (if any)
  also return Sims (if any)
  @param fileList: list of input files
  @param statistic: instance of coincStatistic, to use in creating coincs
  """
  if not fileList:
    return coincInspiralTable(), None

  if not (isinstance(statistic,coincStatistic)):
    raise TypeError, "invalid statistic, must be coincStatistic"

  sims = None
  coincs = None

  for thisFile in fileList:
    doc = utils.load_filename(thisFile, gz = (thisFile or "stdin").endswith(".gz"))
    # extract the sim inspiral table
    try: 
      simInspiralTable = \
          table.get_table(doc, lsctables.SimInspiralTable.tableName)
      if sims: sims.extend(simInspiralTable)
      else: sims = simInspiralTable
    except: simInspiralTable = None

    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: snglInspiralTable = None
    if snglInspiralTable:
      coincFromFile = coincInspiralTable(snglInspiralTable,statistic)
      if simInspiralTable: 
        coincFromFile.add_sim_inspirals(simInspiralTable) 
      if coincs: coincs.extend(coincFromFile)
      else: coincs = coincFromFile
  return coincs, sims


########################################
class coincStatistic:
  """
  This class specifies the statistic to be used when dealing with coincident events.
  It also contains parameter for such cases as the BBH bitten-L statistics.
  """

  __slots__ = ["name","a","b","rsq","bl", "eff_snr_denom_fac"]

  def __init__(self, name, a=0, b=0, denom_fac=250.0):
    self.name=name
    self.a=a
    self.b=b
    self.rsq=0
    self.bl=0
    self.eff_snr_denom_fac = denom_fac

  def __str__(self):
    return self.name

  def get_bittenl(self, bl, snr ):
    blx=self.a*snr-self.b
    if bl==0:
      return blx
    else:
      return min(bl, blx)


#######################################
class coincInspiralTable:
  """
  Table to hold coincident inspiral triggers.  Coincidences are reconstructed 
  by making use of the event_id contained in the sngl_inspiral table.
  The coinc is a dictionary with entries: event_id, numifos, stat, and 
  each available IFO (G1, H1, etc.).
  The stat is set by default to the snrsq: the sum of the squares of the snrs 
  of the individual triggers.
  """
  class row(object):
    __slots__ = ["event_id", "numifos", "stat", "likelihood",
                 "sim", "rsq", "bl", "fap"] + list(ifos)
    
    def __init__(self, event_id, numifos = 0, stat = 0, likelihood = 0):
      self.event_id = event_id
      self.numifos = numifos
      self.stat = stat
      self.likelihood = likelihood
      self.rsq=0
      self.bl=0
      self.fap = 0.0
      
    def add_trig(self,trig,statistic):
      # Coincidence IDs are intended to be unique.  If there is a collision,
      # multiple triggers from the same ifo can get mixed together.  This is
      # a serious problem.  This won't detect all cases, but with more and
      # more triggers being added, it becomes increasingly likely that
      # we'll notice and halt the program.
      assert not hasattr(self, trig.ifo), "Trying to add %s trigger to a"\
        " coincidence for the second time. Coincidence so far:\n%s"\
        "\n\nTrigger:\n%s" % (trig.ifo, dict([(x, getattr(self, x)) for x in \
        self.__slots__ if hasattr(self, x)]), trig.event_id)
      
      self.numifos +=1
      if statistic.name == 'effective_snr':
        self.stat = (self.stat**2 + trig.get_effective_snr(statistic.eff_snr_denom_fac)**2)**(1./2)      
      elif 'bitten_l' in statistic.name:
        snr=trig.snr
        self.rsq= (self.rsq**2 + snr**2)**(1./2)
        self.bl=statistic.get_bittenl( self.bl, snr )
        self.stat=min( self.bl, self.rsq )
        if statistic.name == 'bitten_lsq' and self.numifos >2:
          self.stat = self.rsq
      elif 'far' == statistic.name:
        self.stat = trig.get_far()
      elif 'ifar' == statistic.name:
        self.stat = trig.get_ifar()
      elif 'lvS5stat' == statistic.name:
        self.stat = trig.get_lvS5stat()
      else:
        self.stat = (self.stat**2 + getattr(trig,statistic.name)**2)**(1./2)
      
      # sets the data for the single inspiral trigger
      setattr(self,trig.ifo,trig)
      
    def add_sim(self,sim):
      setattr(self,"sim",sim)

    def get_ifos(self): 
      ifo_string = ""
      ifolist_in_coinc = []
      for ifo in ifos:
        if hasattr(self,ifo):
          ifo_string = ifo_string + ifo
          ifolist_in_coinc.append(ifo)

      return ifo_string, ifolist_in_coinc
    
    def _get_slide_num(self):
      slide_num = (int(self.event_id) % 1000000000) // 100000
      if slide_num > 5000: slide_num = 5000 - slide_num
      return slide_num
    slide_num = property(fget=_get_slide_num)

    def get_time( self ):
      for ifo in ifos:
        if hasattr(self,ifo):
          return getattr(self, ifo).end_time+getattr(self, ifo).end_time_ns*1.0e-9
      raise ValueError, "This coincident trigger does not contain any "\
            "single trigger.  This should never happen."

    def get_gps_times( self ):
      """
      Return a dictionary of the GPS times associated with each
      trigger in the coincidence. The time is stored in LIGOTimeGPS format. 
      """
      gpstimes={}
      for ifo in ifos:
        if hasattr(self,ifo):
          gpstimes[ifo]= LIGOTimeGPS(getattr(self, ifo).end_time, getattr(self, ifo).end_time_ns)
      return gpstimes

    def __iter__(self):
      """
      Return an iterator over the triggers in this coinc.
      """
      for ifo in ifos:
        if hasattr(self, ifo):
          yield getattr(self, ifo)

  def __init__(self, inspTriggers = None, stat = None):
    """
    @param inspTriggers: a metaDataTable containing inspiral triggers 
                         from which to construct coincidences
    @param stat:         an instance of coincStatistic
    """
    self.stat = stat
    self.sngl_table = inspTriggers
    self.sim_table = None
    self.rows = []
    if inspTriggers is None:
      return

    # At present, coincidence is recorded by thinca by over-writing event_ids.
    # The event_ids uniquely label coincidences now, rather than triggers.
    row_dict = {}
    unique_id_list = []
    for trig in inspTriggers:
      event_id = trig.event_id
      if event_id not in row_dict:
        unique_id_list.append(event_id)
        row_dict[event_id] = self.row(event_id)
      row_dict[event_id].add_trig(trig, stat)

    # make sure that there are at least two ifos in each coinc; restore order
    pruned_rows = [row_dict[k] for k in unique_id_list \
      if row_dict[k].numifos > 1]

    self.rows = pruned_rows
    
  def __len__(self):
    return len(self.rows)
  
  def append(self,row):
    self.rows.append(row)

  def extend(self,rows):
    self.rows.extend(rows)

  def __getitem__(self, i):
    """
    Retrieve the value in this column in row i.
    """
    return self.rows[i]

  def getstat(self):
    return numpy.array([c.stat for c in self.rows], dtype=float)

  def sort(self, descending = True):
    """
    Sort the list based on stat value 
    default is to descending
    """
    stat_list = [ (coinc.stat, coinc) for coinc in self.rows ]
    stat_list.sort()
    if descending:
      stat_list.reverse()
    self.rows = [coinc for (stat,coinc) in stat_list]

  def calculate_fap(self, stats, use_likelihood = False):
    """
    Calculates false alarm probability for each coinc using stats array.
    @param stats: array of loudest statistics forund in each of the time slides
    """

    for coinc in self:
      if not use_likelihood:
        index = numpy.searchsorted(stats, coinc.stat)
        if index == 0:
          coinc.fap = 100.0
        elif index == len(stats):
          coinc.fap = 0.0
        else:
          coinc.fap = 100.0 * float((len(stats) - index)) /float(len(stats))

  def getslide(self, slide_num):
    """
    Return the triggers with a specific slide number.
    @param slide_num: the slide number to recover (contained in the event_id)
    """
    slide_coincs = coincInspiralTable(stat=self.stat)
    slide_coincs.sngl_table = self.sngl_table
    slide_coincs.extend([c for c in self.rows if c.slide_num == slide_num])
    return slide_coincs 

  def coincinclude(self, ifolist):
    """
    Return the coincs which have triggers from the ifos in ifolist.
    @param ifolist: a list of ifos 
    """
    selected_coincs = coincInspiralTable(stat=self.stat)
    selected_coincs.sngl_table = self.sngl_table
    for coinc in self:
      keep_trig = True
      for ifo in ifolist:
        if hasattr(coinc,ifo) == False:
          keep_trig = False
          break
            
      if keep_trig == True:
        selected_coincs.append(coinc)
        
    return selected_coincs

  def coinctype(self, ifolist):
    """
    Return the coincs which are from ifos.
    @param ifolist: a list of ifos 
    """
    coincs = self.coincinclude(ifolist)
    selected_coincs = coincInspiralTable(stat=self.stat)
    selected_coincs.sngl_table = self.sngl_table
    for coinc in coincs:
      if coinc.numifos == len(ifolist):
        selected_coincs.append(coinc)
        
    return selected_coincs

  def removecoinctype(self, ifolist):
    """
    Return the coincs which are NOT from the coinc type made from combining all
    the ifos in the ifolist.
    @param ifolist: a list of ifos 
    """
    ifolist.sort()
    selected_coincs = coincInspiralTable(stat=self.stat)
    selected_coincs.sngl_table = self.sngl_table
    for coinc in self:
      thiscoinctype = coinc.get_ifos()[1]
      thiscoinctype.sort()
      if not (ifolist == thiscoinctype):
        selected_coincs.append(coinc)

    return selected_coincs

  def getsngls(self, ifo):
    """
    Return the sngls for a specific ifo.
    @param ifo: ifo for which to retrieve the single inspirals.
    """
    from glue.ligolw import table 
    try: ifoTrigs = table.new_from_template(self.sngl_table)
    except: ifoTrigs = lsctables.New(lsctables.SnglInspiralTable)
    for coinc in self:
      if hasattr(coinc,ifo): 
        ifoTrigs.append(getattr(coinc,ifo))
        
    return ifoTrigs

  def vetoed(self, seglist):
    """
    Returns a list of coincident triggers that are vetoed
    entirely by a segment of the segment list.
    A coincident trigger is added to the list only
    if all triggers lie within a segment. If any
    trigger lies outside any segment, it is not added.
    @param seglist: segment list used to veto coincidences
    """
    vetoed_coincs = coincInspiralTable(stat=self.stat)

    # loop over all the coincident triggers
    for coinc in self:
      flagVeto = True
      for ifo in ifos:
        if hasattr(coinc, ifo):
          # if any trigger is outside, do not add
          if getattr(coinc,ifo).get_end() not in seglist:
            flagVeto = False
            break

      # append to list of vetoed coincs if and only if
      # all triggers lie within the seglist
      if flagVeto:
        vetoed_coincs.append( coinc )

    return vetoed_coincs

  def cluster_core(self, cluster_window):
    """
    Return the clustered triggers, returning the one with the largest stat in 
    each fixed cluster_window, exactly as in lalapps_coire
    (in package/tools/CoincInspiralUtils.c:XLALClusterCoincInspiralTable)
    
    @param cluster_window: length of time over which to cluster (seconds)
    """

    # first we need to time-sort the coincidences
    stat_list = [ (coinc.get_time(), coinc) for coinc in self.rows ]
    stat_list.sort()
    self.rows = [coinc for (t,coinc) in stat_list]

    # initialize some indices (could work with just one)
    # but with two its easier to read
    thisCoinc = 0
    nextCoinc = 1
    while nextCoinc<len(self.rows):

      # get the time for both indices
      thisTime = self.rows[thisCoinc].get_time()
      nextTime = self.rows[nextCoinc].get_time()

      # are the two coincs within the time-window?
      if nextTime-thisTime<cluster_window:

        # get the statistic values
        thisStat = self.rows[thisCoinc].stat
        nextStat = self.rows[nextCoinc].stat

        # and remove the coinc which has the lower coinc value
        if (nextStat>thisStat):          
          self.rows.pop(thisCoinc)      
        else:
          self.rows.pop(nextCoinc)

      else:
        # the two coincidences are NOT in the time-window
        # so must increase index
        thisCoinc+=1
        nextCoinc+=1
        

  def cluster(self, cluster_window, numSlides = None):
    """
    New clustering method working the same way as the clustering method
    used in lalapps_coire
    (in package/tools/CoincInspiralUtils.c:XLALClusterCoincInspiralTable)

    @param cluster_window: length of time over which to cluster (seconds)
    @param numSlides: number of slides if this are time-slide triggers
    """
    
    if not numSlides:
      
      # just cluster them
      self.cluster_core(cluster_window)

    else:
      # create a new coincInspiralTable
      cluster = coincInspiralTable()

      # need to cluster all triggers from each slide individually
      for slide in range(-numSlides, numSlides+1):

        # get the coincs belonging to slide 'slide'
        slideCoinc = self.getslide( slide )

        # and cluster these
        slideCoinc.cluster_core( cluster_window )

        # add the clustered coincs to 'cluster'
        cluster.extend( slideCoinc )

      # replace the clustered triggers
      self.rows = cluster.rows

    # return this object itself
    return self
      
  
  def add_sim_inspirals(self,sim_inspiral):
    """
    FIXME: We should really store the sim coincidence info in the event_id
    Method to add simulated injections to a list of coincs

    @param sim_inspiral: a simInspiralTable
    """
    self.sim_table = sim_inspiral
    # check that the number of sims matches the number of coincs:
    if len(self) != len(sim_inspiral):
      raise ValueError, "Number of injections doesn't match number of coincs"

    for i in range(len(self)):
      self[i].add_sim(sim_inspiral[i])


  def add_missed_sims(self,sim_inspiral):
    """
    Add missed sim inspirals to the list of coincs, set the stat = -1
    @param sim_inspiral: a simInspiralTable
    """
    for sim in sim_inspiral:
      row = coincInspiralTable.row(-1,stat=-1)
      row.add_sim(sim)
      self.append(row)

  def return_sim_inspirals(self,statistic = None,thresh = 0):
    """
    Method to return the sim_inspiral table associated to the coincs.
    If thresh is specified, only return sims from those coincs whose stat
    exceeds thresh (or is under thresh if statistic == far).

    @param statistic: the statistic to use
    @param thresh: the threshold on the statistic
    """
    from glue.ligolw import table 
    try: simInspirals = table.new_from_template(self.sim_table)
    except: simInspirals = lsctables.New(lsctables.SimInspiralTable)
    for coinc in self:
      if statistic == 'far':
        if (hasattr(coinc,"sim")) and \
            (((coinc.stat <= thresh) and (coinc.stat >= 0)) or (thresh < 0)):
          simInspirals.append(coinc.sim)
      else:
        if (hasattr(coinc,"sim")) and (coinc.stat >= thresh):
          simInspirals.append(coinc.sim)
    
    return simInspirals
    
  def partition_by_stat(self, threshold):
    """
    Return (triggers with stat < threshold,
    triggers with stat == threshold,
    triggers with stat > threshold).

    The set of (stat == threshold) is of zero measure, but often, as when
    doing something like the loudest event statistic, the threshold is taken
    from a coinc in self.
    """
    stats = self.getstat()

    lesser_coincs = coincInspiralTable(stat=self.stat)
    lesser_coincs.extend([self[i] for i in (stats < threshold).nonzero()[0]])

    equal_coincs = coincInspiralTable(stat=self.stat)
    equal_coincs.extend([self[i] for i in (stats == threshold).nonzero()[0]])

    greater_coincs = coincInspiralTable(stat=self.stat)
    greater_coincs.extend([self[i] for i in (stats > threshold).nonzero()[0]])

    return lesser_coincs, equal_coincs, greater_coincs

  def getTotalMass(self,mLow,mHigh):
    """
    Return triggers with mLow <= mean total mass < mHigh
    @param mLow: a float
    @param mHigh: a float
    """
    triggers_in_mass_range = coincInspiralTable()
    for coinc in self:
      mass_numer = 0.0
      mass_denom = float(coinc.numifos)
      for ifo in ifos:
        if hasattr(coinc,ifo):
          mass_numer += getattr(coinc,ifo).mass1 + getattr(coinc,ifo).mass2
      mean_total_mass = mass_numer / mass_denom
      if (mean_total_mass >= mLow) and (mean_total_mass < mHigh):
        triggers_in_mass_range.append(coinc)

    return triggers_in_mass_range

  def getChirpMass(self,mLow,mHigh):
    """
    Return triggers with mLow <= mean chirp mass < mHigh
    @param mLow: a float
    @param mHigh: a float
    """
    triggers_in_mass_range = coincInspiralTable()
    for coinc in self:
      mass_numer = 0.0
      mass_denom = float(coinc.numifos)
      for ifo in ifos:
        if hasattr(coinc,ifo):
          mass_numer += getattr(coinc,ifo).mchirp
      mean_total_mass = mass_numer / mass_denom
      if (mean_total_mass >= mLow) and (mean_total_mass < mHigh):
        triggers_in_mass_range.append(coinc)

    return triggers_in_mass_range

  def getEThincaValues(self,ifos):
    """
    Return ethinca values for the coincidences
    @param ifos: a list of the 2 ifos
    """
    ethinca = numpy.zeros(len(self), dtype=float)
    for i,coinc in enumerate(self):
      if hasattr(coinc, ifos[0]) and hasattr(coinc, ifos[1]):
        try: 
          ethinca[i] = XLALCalculateEThincaParameter(getattr(coinc, ifos[0]),
                                                     getattr(coinc, ifos[1]))
        except ValueError:
          ethinca[i] = 2.0
          
    return ethinca

  def getSimpleEThincaValues(self, ifos):
    """
    Return simple ethinca values for the coincidences of the ifo pair ifos.
    
    For coincs that do not have both ifos specified, return 0.0.
    """
    ethinca = numpy.zeros(len(self), dtype=float)
    for i,coinc in enumerate(self):
      if hasattr(coinc, ifos[0]) and hasattr(coinc, ifos[1]):
        ethinca[i] = simpleEThinca(getattr(coinc, ifos[0]),
                                   getattr(coinc, ifos[1]))
    return ethinca

	
  def getTriggersWithinEpsilon(self, candidate, epsilon, epsilon_ball_type = "version2"):
	"""
	Return distance squared between candidate and coincident triggers
	using the following metric
	d^2 = ( 1/n ) * sum_i ( |cparam_i - trigparam_i|^2 / cparam_i^2 )
	@param candidate: a coincInspiral describing a candidate
	"""
	# epsilon ball version 1
	##############################################################
	if epsilon_ball_type == "version1":
	  triggers_within_epsilon = coincInspiralTable()
	  epsilon_sq = epsilon * epsilon

	  # setting which dimension should be used
	  dim_stat = False
	  dim_mchirp_err = False
	  dim_average_ethinca = False
	  dim_effective_distance = False
	  
	  # setting scales for dimensions
	  stat_scale = 1.0
	  mchirp_err_scale = 0.1 
	  average_ethinca_scale = 0.1
	  eff_dist_diff_scale = 0.1
	  
	  # calculating parameters of the candidate 
	  c_ifos,ifolist = candidate.get_ifos()
	  # combined effective snr
	  if dim_stat:
		c_stat = candidate.stat
	  # combined fractional diff in chirp mass
	  if dim_mchirp_err:
		c_mchirp_err_sq = 0.0
	  # combined ethinca
	  if dim_average_ethinca:
		c_sum_ethinca = 0.0
	  # combined fractional diff in effective distance
	  if dim_effective_distance:
	   c_eff_dist_diff_sq = 0.0


	  counter = 0.0
	  for ifo1 in ifolist:
		for ifo2 in ifolist:
		  if ifo1 < ifo2:
			counter += 1.0
			c_ifo1 = getattr(candidate, ifo1)
			c_ifo2 = getattr(candidate, ifo2)
			if dim_mchirp_err:
			  c_mchirp_err_sq += (2.0*(c_ifo2.mchirp - c_ifo1.mchirp)/(c_ifo2.mchirp + c_ifo1.mchirp))**2
			if dim_average_ethinca:
			  c_sum_ethinca += simpleEThinca(c_ifo1, c_ifo2)
			if dim_effective_distance:
			  c_eff_dist_diff_sq += (2.0*(c_ifo2.eff_distance - c_ifo1.eff_distance)/(c_ifo2.eff_distance + c_ifo1.eff_distance))**2
	  if dim_mchirp_err:    
		c_mchirp_err = c_mchirp_err_sq**(0.5)
	  if dim_average_ethinca:
		c_average_ethinca = c_sum_ethinca/counter
	  if dim_effective_distance:
		c_eff_dist_diff = c_eff_dist_diff_sq**(0.5)

	  # loop over triggers
	  for trig in self:
		trig_ifos,tmplist = trig.get_ifos()
		if c_ifos == trig_ifos:
		  tmp_d_squared = 0.0

		  # distance^2 apart in combined effective snr
		  if dim_stat:
			t_stat = trig.stat
			tmp_d_squared += ((c_stat - t_stat)/stat_scale)**2

		  if dim_mchirp_err:
			t_mchirp_err_sq = 0.0
		  if dim_average_ethinca:
			t_sum_ethinca = 0.0
		  if dim_effective_distance:
			t_eff_dist_diff_sq = 0.0
	
		  counter = 0.0
		  for ifo1 in ifolist:
			for ifo2 in ifolist:
			  if ifo1 < ifo2:
				counter += 1.0
				t_ifo1 = getattr(trig, ifo1)
				t_ifo2 = getattr(trig, ifo2)
				if dim_mchirp_err:
				  t_mchirp_err_sq += (2.0*(t_ifo2.mchirp - t_ifo1.mchirp)/(t_ifo2.mchirp + t_ifo1.mchirp))**2
				if dim_average_ethinca:
				  t_sum_ethinca += simpleEThinca(t_ifo1, t_ifo2)
				if dim_effective_distance:
				  t_eff_dist_diff_sq += (2.0*(t_ifo2.eff_distance - t_ifo1.eff_distance)/(t_ifo2.eff_distance + t_ifo1.eff_distance))**2

		  # distance^2 apart in combined fractional difference in mchirp
		  if dim_mchirp_err:           
			t_mchirp_err = t_mchirp_err_sq**(0.5)
			tmp_d_squared += ((c_mchirp_err - t_mchirp_err)/mchirp_err_scale)**2
		  # distance^2 apart in combined ethinca
		  if dim_average_ethinca:
			t_average_ethinca = t_sum_ethinca/counter
			tmp_d_squared += ((c_average_ethinca - t_average_ethinca)/average_ethinca_scale)**2

		  # distance^ apart in combined fractional difference in effective distance
		  if dim_effective_distance:
			t_eff_dist_diff = t_eff_dist_diff_sq**(0.5)
			tmp_d_squared += ((c_eff_dist_diff - t_eff_dist_diff)/eff_dist_diff_scale)**2
		  if ( tmp_d_squared < epsilon_sq ):
			triggers_within_epsilon.append(trig)
			
	# epsilon ball version2 
	#######################################################################
	if epsilon_ball_type == "version2":
	  triggers_within_epsilon = coincInspiralTable()
	  epsilon_sq = epsilon * epsilon
	  # setting which dimension should be used
	  dim_eff_snr = True
	  dim_mchirp = True
	  dim_ethinca = False
	  dim_eff_distance = True
	  
	  # setting scales for dimensions
	  #eff_snr_scale = 1.0
	  #mchirp_scale = 0.1 
	  ethinca_scale = 0.1
	  eff_dist_theta_scale = 2.0*cmath.pi/float(180)
	  eff_dist_phi_scale =  10.0*cmath.pi/float(180)
	  mchirp_theta_scale = 2.0*cmath.pi/float(180)
	  mchirp_phi_scale = 180.0*cmath.pi/float(180) 
	
	  # calculating parameters of the candidate 
	  c_ifos,ifolist = candidate.get_ifos()
	  
	  # if candidate is  a double 
	  if len(ifolist) == 2:
		c_ifo1 = getattr(candidate, ifolist[0])
		c_ifo2 = getattr(candidate, ifolist[1])
		if dim_eff_snr:
		  c_eff_snr_ifo1 = c_ifo1.get_effective_snr()
		  c_eff_snr_ifo2 = c_ifo2.get_effective_snr()
		if dim_eff_distance:
		  c_D_eff_rms,c_eff_dist_theta = convert_to_polar_coord(c_ifo1.eff_distance, c_ifo2.eff_distance)
		if dim_mchirp:
		  c_mchirp_rms, c_mchirp_theta = convert_to_polar_coord(c_ifo1.mchirp, c_ifo2.mchirp)
		if dim_ethinca:
		  c_ethinca = simpleEThinca(c_ifo1, c_ifo2)
                print "eff_snr_dim:", str(c_eff_snr_ifo1), str(c_eff_snr_ifo2)
                print "D_eff_dim:", str(c_D_eff_rms), str(c_eff_dist_theta)
                print "mchirp_dim:", str(c_mchirp_rms), str(c_mchirp_theta)
	  # if candidate is a triple	  
	  if len(ifolist) == 3:
		c_ifo1 = getattr(candidate, ifolist[0])
		c_ifo2 = getattr(candidate, ifolist[1])
		c_ifo3 = getattr(candidate, ifolist[2])
		if dim_eff_snr:
		  c_eff_snr_ifo1 = c_ifo1.get_effective_snr()
		  c_eff_snr_ifo2 = c_ifo2.get_effective_snr()
		  c_eff_snr_ifo3 = c_ifo3.get_effective_snr()
		if dim_eff_distance:
		  c_D_eff_rms, c_eff_dist_theta, c_eff_dist_phi = convert_to_spherical_coord(c_ifo1.eff_distance, c_ifo2.eff_distance, c_ifo3.eff_distance)
		if dim_mchirp:
		  c_mchirp_rms, c_mchirp_theta, c_mchirp_phi = convert_to_spherical_coord(c_ifo1.mchirp, c_ifo2.mchirp, c_ifo3.mchirp)
		if dim_ethinca:
		  c_ethinca_12 = simpleEThinca(c_ifo1, c_ifo2)
		  c_ethinca_13 = simpleEThinca(c_ifo1, c_ifo3)
		  c_ethinca_23 = simpleEThinca(c_ifo2, c_ifo3)
		  
	  # loop over triggers
	  for trig in self:
		trig_ifos,tmplist = trig.get_ifos()
		if c_ifos == trig_ifos:
		# if candidate is  a double 
		  if len(ifolist) == 2:
			t_ifo1 = getattr(trig, ifolist[0])
			t_ifo2 = getattr(trig, ifolist[1])
			if dim_eff_snr:
			  t_eff_snr_ifo1 = t_ifo1.get_effective_snr()
			  t_eff_snr_ifo2 = t_ifo2.get_effective_snr()
			if dim_eff_distance:
			  t_D_eff_rms,t_eff_dist_theta = convert_to_polar_coord(t_ifo1.eff_distance, t_ifo2.eff_distance)
			if dim_mchirp:
			  t_mchirp_rms, t_mchirp_theta = convert_to_polar_coord(t_ifo1.mchirp, t_ifo2.mchirp)
			if dim_ethinca:
			  t_ethinca = simpleEThinca(t_ifo1, t_ifo2)
			score = 0.0
			if not dim_eff_snr or ((frac_error_sq(c_eff_snr_ifo1, t_eff_snr_ifo1) < epsilon_sq) and (frac_error_sq(c_eff_snr_ifo2, t_eff_snr_ifo2) < epsilon_sq)):
			  score += 1.0	
			if not dim_eff_distance or ((frac_error_sq(c_D_eff_rms, t_D_eff_rms) < epsilon_sq) and (abs(c_eff_dist_theta - t_eff_dist_theta) < eff_dist_theta_scale)):
			  score += 1.0
                        #WARNING: epsilon for mchirp is hardcoded to be 0.5 
			if not dim_mchirp or ((frac_error_sq(c_mchirp_rms, t_mchirp_rms) < 0.25) and (abs(c_mchirp_theta - t_mchirp_theta) < mchirp_theta_scale)):
			  score += 1.0
			if not dim_ethinca or (abs(c_ethinca - t_ethinca) < ethinca_scale):
			  score += 1.0
			if score == 4.0:
			  triggers_within_epsilon.append(trig)
                          print "candidate:" , str(c_ifo1.event_id)
                          print "eff_snr_dim:", str(c_eff_snr_ifo1), str(c_eff_snr_ifo2)
                          print "D_eff_dim:", str(c_D_eff_rms), str(c_eff_dist_theta)
                          print "mchirp_dim:", str(c_mchirp_rms), str(c_mchirp_theta)
                          print "score=", score
                          print "trigger:", str(t_ifo1.event_id)
                          print "eff_snr_dim:", str(t_eff_snr_ifo1), str(t_eff_snr_ifo2)
                          print "D_eff_dim:", str(t_D_eff_rms), str(t_eff_dist_theta)
                          print "mchirp_dim:", str(t_mchirp_rms), str(t_mchirp_theta)
                          print "epsilon_sq=", epsilon_sq
                          print "snr condition:", str(frac_error_sq(c_eff_snr_ifo1, t_eff_snr_ifo1)), str(frac_error_sq(c_eff_snr_ifo2, t_eff_snr_ifo2))
                          print "eff distance condition:", str(frac_error_sq(c_D_eff_rms, t_D_eff_rms)), str(abs(c_eff_dist_theta - t_eff_dist_theta)), str(eff_dist_theta_scale) 
                          print "mchirp:", str(frac_error_sq(c_mchirp_rms, t_mchirp_rms)), str(abs(c_mchirp_theta - t_mchirp_theta)), str(mchirp_theta_scale) 
			  
		  # if candidate is a triple	  
		  if len(ifolist) == 3:
			t_ifo1 = getattr(trig, ifolist[0])
			t_ifo2 = getattr(trig, ifolist[1])
			t_ifo3 = getattr(trig, ifolist[2])
			if dim_eff_snr:
			  t_eff_snr_ifo1 = t_ifo1.get_effective_snr()
			  t_eff_snr_ifo2 = t_ifo2.get_effective_snr()
			  t_eff_snr_ifo3 = t_ifo3.get_effective_snr()
			if dim_eff_distance:
			  t_D_eff_rms, t_eff_dist_theta, t_eff_dist_phi = convert_to_spherical_coord(t_ifo1.eff_distance, t_ifo2.eff_distance, t_ifo3.eff_distance)
			if dim_mchirp:
			  t_mchirp_rms, t_mchirp_theta, t_mchirp_phi = convert_to_spherical_coord(t_ifo1.mchirp, t_ifo2.mchirp, t_ifo3.mchirp)
			if dim_ethinca:
			  t_ethinca_12 = simpleEThinca(t_ifo1, t_ifo2)
			  t_ethinca_13 = simpleEThinca(t_ifo1, t_ifo3)
			  t_ethinca_23 = simpleEThinca(t_ifo2, t_ifo3)
			score = 0.0
			if not dim_eff_snr or ((frac_error_sq(c_eff_snr_ifo1, t_eff_snr_ifo1) < epsilon_sq) and (frac_error_sq(c_eff_snr_ifo2, t_eff_snr_ifo2) < epsilon_sq) and (frac_error_sq(c_eff_snr_ifo3, t_eff_snr_ifo3) < epsilon_sq)):
			  score += 1.0	
			if not dim_eff_distance or ((frac_error_sq(c_D_eff_rms, t_D_eff_rms) < epsilon_sq) and (abs(c_eff_dist_theta - t_eff_dist_theta) < eff_dist_theta_scale) and (abs(c_eff_dist_phi - t_eff_dist_phi) < eff_dist_phi_scale)):
			  score += 1.0
                        #WARNING: epsilon for mchirp is hardcoded to be 0.5
			if not dim_mchirp or ((frac_error_sq(c_mchirp_rms, t_mchirp_rms) < 0.25) and (abs(c_mchirp_theta - t_mchirp_theta) < mchirp_theta_scale) and (abs(c_mchirp_phi - t_mchirp_phi) < mchirp_phi_scale)):
			  score += 1.0
			if not dim_ethinca or ((abs(c_ethinca_12 - t_ethinca_12) < ethinca_scale) and (abs(c_ethinca_13 - t_ethinca_13) < ethinca_scale) and (abs(c_ethinca_23 - t_ethinca_23) < ethinca_scale)):
			  score += 1.0
			if score == 4.0:
			  triggers_within_epsilon.append(trig)

		
		
		


	if epsilon_ball_type == "IFAR":
	  triggers_within_epsilon = coincInspiralTable()
	  
	  # calculating bins for mchirp
	  mchirp1 = 2.0*(0.5)**(6.0/5.0)
	  mchirp2 = 8.0*(0.5)**(6.0/5.0)
	  mchirp3 = 17.0*(0.5)**(6.0/5.0)
	  mchirp4 = 35.0*(0.5)**(6.0/5.0)
	  
	  # calculating parameters of the candidate 
	  c_ifos,ifolist = candidate.get_ifos()
	  c_mchirp_total = 0.0
	  
	  counter_sngl = 0.0
	  for ifo in ifolist:
		c_sngl = getattr(candidate, ifo) 
		c_mchirp_total += c_sngl.mchirp
		counter_sngl += 1.0
	  c_mchirp_average = c_mchirp_total/counter_sngl
	  c_eff_snr = candidate.stat

	  # loop over triggers
	  for trig in self:
		trig_ifos,tmplist = trig.get_ifos()
		if c_ifos == trig_ifos:
		  t_mchirp_total = 0.0
		  counter_sngl = 0.0  
		  for ifo in ifolist:
			t_sngl = getattr(trig, ifo)
			t_mchirp_total += t_sngl.mchirp
			counter_sngl += 1.0
		  t_mchirp_average = t_mchirp_total/counter_sngl
		  t_eff_snr = trig.stat
		  if t_eff_snr >= c_eff_snr:
			if (mchirp1 <= c_mchirp_average < mchirp2) and (mchirp1 <= t_mchirp_average < mchirp2):
			  triggers_within_epsilon.append(trig)
			elif (mchirp2 <= c_mchirp_average < mchirp3) and (mchirp2 <= t_mchirp_average < mchirp3):
			  triggers_within_epsilon.append(trig)
			elif (mchirp3 <= c_mchirp_average < mchirp4) and (mchirp3 <= t_mchirp_average <= mchirp4): 
			  triggers_within_epsilon.append(trig)


        #for ifo1 in ifolist: 
          #c = getattr(candidate,ifo1)
          #t = getattr(trig,ifo1)
          
          # distance ^2 apart in combined effective snr (units 1 snr)
          #c_lambda = c.get_effective_snr()
          #t_lambda = t.get_effective_snr()
          #tmp_d_squared += (c_lambda - t_lambda)**2
          #param_counter += 1           
          # distance^2 apart in snr
          #c_lambda = c.snr
          #t_lambda = t.snr
          #tmp_d_squared += ( 1.0 - t_lambda / c_lambda )**2
          #param_counter += 1

          # distance^2 apart in chi squared
          #c_lambda = c.chisq
          #t_lambda = t.chisq
          #tmp_d_squared += ( 1.0 - t_lambda / c_lambda )**2
          #param_counter += 1

          # distance^2 apart in mchirp
          #c_lambda = c.mchirp
          #t_lambda = t.mchirp
          #tmp_d_squared += ( 1.0 - t_lambda / c_lambda )**2
          #param_counter += 1

          # distance^2 apart in effective distance
          #c_lambda = c.eff_distance
          #t_lambda = t.eff_distance
          #tmp_d_squared += ( 1.0 - t_lambda / c_lambda )**2
          #param_counter += 1

          # distance^2 apart in ethinca
          #for ifo2 in ifolist:
          #  if ifo1 < ifo2:
          #    c_lambda = simpleEThinca(c, getattr(candidate, ifo2))
          #    t_lambda = simpleEThinca(t, getattr(trig, ifo2))
              
              # check if the ethinca parameter of the candidate is not zero
          #    if c_lambda > 0.0:
          #      sEThinca_scale = c_lambda
          #    else:
          #      sEThinca_scale = 0.5
              #tmp_d_squared += ( 1.0 - t_lambda / c_lambda )**2
              # FIX ME! We need to make a choice for the sEThinca_scale that sets the scale for SimpleEThinca parameter 
          #    tmp_d_squared += ((c_lambda - t_lambda) / sEThinca_scale)**2
          #    param_counter += 1
        #if ( tmp_d_squared < epsilon_sq ):
        #  triggers_within_epsilon.append(trig)

	return triggers_within_epsilon
  
  def getTriggersWithinEpsilonBox(self, candidate, dim_eff_snr, dim_mchirp, dim_ethinca, dim_eff_distance, eff_snr_epsilon, eff_dist_epsilon, mchirp_epsilon,\
	eff_dist_theta_scale, eff_dist_phi_scale, mchirp_theta_scale, mchirp_phi_scale, ethinca_scale):
	""" 
	Returns triggers that lie in the parameter region centered around the candidate. 
	Region is defined by the set of conditions. 
	"""
	triggers_within_epsilon = coincInspiralTable()
	eff_snr_epsilon_sq = eff_snr_epsilon * eff_snr_epsilon
	eff_dist_epsilon_sq = eff_dist_epsilon * eff_dist_epsilon
	mchirp_epsilon_sq = mchirp_epsilon * mchirp_epsilon
	# setting which dimension should be used
	#dim_eff_snr = True
	#dim_mchirp = True
	#dim_ethinca = False
	#dim_eff_distance = True
	
	# setting scales for dimensions
	#eff_snr_scale = 1.0
	#mchirp_scale = 0.1 
	#ethinca_scale = 0.1
	#eff_dist_theta_scale = 2.0*cmath.pi/float(180)
	#eff_dist_phi_scale =  10.0*cmath.pi/float(180)
	#mchirp_theta_scale = 2.0*cmath.pi/float(180)
	#mchirp_phi_scale = 180.0*cmath.pi/float(180) 
  
	# calculating parameters of the candidate 
	c_ifos,ifolist = candidate.get_ifos()
	
	# if candidate is  a double 
	if len(ifolist) == 2:
	  c_ifo1 = getattr(candidate, ifolist[0])
	  c_ifo2 = getattr(candidate, ifolist[1])
	  if dim_eff_snr:
		c_eff_snr_ifo1 = c_ifo1.get_effective_snr()
		c_eff_snr_ifo2 = c_ifo2.get_effective_snr()
	  if dim_eff_distance:
		c_D_eff_rms,c_eff_dist_theta = convert_to_polar_coord(c_ifo1.eff_distance, c_ifo2.eff_distance)
	  if dim_mchirp:
		c_mchirp_rms, c_mchirp_theta = convert_to_polar_coord(c_ifo1.mchirp, c_ifo2.mchirp)
	  if dim_ethinca:
		c_ethinca = simpleEThinca(c_ifo1, c_ifo2)
	# if candidate is a triple	  
	if len(ifolist) == 3:
	  c_ifo1 = getattr(candidate, ifolist[0])
	  c_ifo2 = getattr(candidate, ifolist[1])
	  c_ifo3 = getattr(candidate, ifolist[2])
	  if dim_eff_snr:
		c_eff_snr_ifo1 = c_ifo1.get_effective_snr()
		c_eff_snr_ifo2 = c_ifo2.get_effective_snr()
		c_eff_snr_ifo3 = c_ifo3.get_effective_snr()
	  if dim_eff_distance:
		c_D_eff_rms, c_eff_dist_theta, c_eff_dist_phi = convert_to_spherical_coord(c_ifo1.eff_distance, c_ifo2.eff_distance, c_ifo3.eff_distance)
	  if dim_mchirp:
		c_mchirp_rms, c_mchirp_theta, c_mchirp_phi = convert_to_spherical_coord(c_ifo1.mchirp, c_ifo2.mchirp, c_ifo3.mchirp)
	  if dim_ethinca:
		c_ethinca_12 = simpleEThinca(c_ifo1, c_ifo2)
		c_ethinca_13 = simpleEThinca(c_ifo1, c_ifo3)
		c_ethinca_23 = simpleEThinca(c_ifo2, c_ifo3)
		
	# loop over triggers
	for trig in self:
	  trig_ifos,tmplist = trig.get_ifos()
	  if c_ifos == trig_ifos:
	  # if candidate is  a double 
		if len(ifolist) == 2:
		  t_ifo1 = getattr(trig, ifolist[0])
		  t_ifo2 = getattr(trig, ifolist[1])
		  if dim_eff_snr:
			t_eff_snr_ifo1 = t_ifo1.get_effective_snr()
			t_eff_snr_ifo2 = t_ifo2.get_effective_snr()
		  if dim_eff_distance:
			t_D_eff_rms,t_eff_dist_theta = convert_to_polar_coord(t_ifo1.eff_distance, t_ifo2.eff_distance)
		  if dim_mchirp:
			t_mchirp_rms, t_mchirp_theta = convert_to_polar_coord(t_ifo1.mchirp, t_ifo2.mchirp)
		  if dim_ethinca:
			t_ethinca = simpleEThinca(t_ifo1, t_ifo2)
		  score = 0.0
		  if not dim_eff_snr or ((frac_error_sq(c_eff_snr_ifo1, t_eff_snr_ifo1) < eff_snr_epsilon_sq) and (frac_error_sq(c_eff_snr_ifo2, t_eff_snr_ifo2) < eff_snr_epsilon_sq)):
			score += 1.0	
		  if not dim_eff_distance or ((frac_error_sq(c_D_eff_rms, t_D_eff_rms) < eff_dist_epsilon_sq) and (abs(c_eff_dist_theta - t_eff_dist_theta) < eff_dist_theta_scale)):
			score += 1.0
		  if not dim_mchirp or ((frac_error_sq(c_mchirp_rms, t_mchirp_rms) < mchirp_epsilon_sq) and (abs(c_mchirp_theta - t_mchirp_theta) < mchirp_theta_scale)):
			score += 1.0
		  if not dim_ethinca or (abs(c_ethinca - t_ethinca) < ethinca_scale):
			score += 1.0
		  if score == 4.0:
			triggers_within_epsilon.append(trig)
			
		# if candidate is a triple	  
		if len(ifolist) == 3:
		  t_ifo1 = getattr(trig, ifolist[0])
		  t_ifo2 = getattr(trig, ifolist[1])
		  t_ifo3 = getattr(trig, ifolist[2])
		  if dim_eff_snr:
			t_eff_snr_ifo1 = t_ifo1.get_effective_snr()
			t_eff_snr_ifo2 = t_ifo2.get_effective_snr()
			t_eff_snr_ifo3 = t_ifo3.get_effective_snr()
		  if dim_eff_distance:
			t_D_eff_rms, t_eff_dist_theta, t_eff_dist_phi = convert_to_spherical_coord(t_ifo1.eff_distance, t_ifo2.eff_distance, t_ifo3.eff_distance)
		  if dim_mchirp:
			t_mchirp_rms, t_mchirp_theta, t_mchirp_phi = convert_to_spherical_coord(t_ifo1.mchirp, t_ifo2.mchirp, t_ifo3.mchirp)
		  if dim_ethinca:
			t_ethinca_12 = simpleEThinca(t_ifo1, t_ifo2)
			t_ethinca_13 = simpleEThinca(t_ifo1, t_ifo3)
			t_ethinca_23 = simpleEThinca(t_ifo2, t_ifo3)
		  score = 0.0
		  if not dim_eff_snr or ((frac_error_sq(c_eff_snr_ifo1, t_eff_snr_ifo1) < eff_snr_epsilon_sq) and (frac_error_sq(c_eff_snr_ifo2, t_eff_snr_ifo2) < eff_snr_epsilon_sq) and (frac_error_sq(c_eff_snr_ifo3, t_eff_snr_ifo3) < eff_snr_epsilon_sq)):
			score += 1.0	
		  if not dim_eff_distance or ((frac_error_sq(c_D_eff_rms, t_D_eff_rms) < eff_dist_epsilon_sq) and (abs(c_eff_dist_theta - t_eff_dist_theta) < eff_dist_theta_scale) and (abs(c_eff_dist_phi - t_eff_dist_phi) < eff_dist_phi_scale)):
			score += 1.0
		  if not dim_mchirp or ((frac_error_sq(c_mchirp_rms, t_mchirp_rms) < mchirp_epsilon_sq) and (abs(c_mchirp_theta - t_mchirp_theta) < mchirp_theta_scale) and (abs(c_mchirp_phi - t_mchirp_phi) < mchirp_phi_scale)):
			score += 1.0
		  if not dim_ethinca or ((abs(c_ethinca_12 - t_ethinca_12) < ethinca_scale) and (abs(c_ethinca_13 - t_ethinca_13) < ethinca_scale) and (abs(c_ethinca_23 - t_ethinca_23) < ethinca_scale)):
			score += 1.0
		  if score == 4.0:
			triggers_within_epsilon.append(trig)


  	return triggers_within_epsilon


		


  def getTriggersInSegment(self, segment):
    """
    Return a new coincInspiralTable with triggers whose end_times lie within
    segment; always use alphabetically first ifo's end_time.
    """
    triggers_within_segment = coincInspiralTable(stat=self.stat)

    for trig in self:
      end_time = getattr(trig, trig.get_ifos()[1][0]).end_time
      if end_time in segment:
        triggers_within_segment.append(trig)

    return triggers_within_segment
	
