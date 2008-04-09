import os
import sys

from numpy import *
from pylab import *

from glue.ligolw import ligolw
from glue.ligolw import utils
from glue.ligolw import table
from glue.ligolw import lsctables
from glue import segments
from glue import segmentsUtils
from pylal import CoincInspiralUtils
from pylal import SimInspiralUtils
from pylal import SnglInspiralUtils
from pylal import tools



#######################################################
## downloadLust
#######################################################

def downloadList( output=None ):
  """
  downloads the latest list of GRB's from
  http://www.uoregon.edu/~ileonor/ligo/s5/grb/online/currentgrbs_html_S5.txt

  Usage:
  downloadList( [filename] )

  If no filename is specified, the output will be placed into
  'recentS5.list'
  """
  if not output:
    output = "recentS5.list"
  
  fetchURL = "http://www.uoregon.edu/~ileonor/ligo/s5/grb/online/"\
             "currentgrbs_html_S5.txt"
  command  = "wget -O %s %s" % (output, fetchURL)
  os.system( command )
  


#######################################################
## massesToTau
#######################################################

def convertMassesToTau( flower, m1, m2 ):
  """
  Convert the given masses to the appropriate values of tau0 and tau3.
  Algorithm stolen from LALInspiralParameterCalc in inspiral/LALInspiralParameterCalc.c
  """
  piFl 	= pi *flower
  MTSUN=4.925e-6
  
  eta = m1*m2/power(m1+m2,2)
  if eta>0.25:
    eta-=1.0e-10
    
  totalMass= (m1+m2) * MTSUN
  
  tau0 	= 5.0/(256.0*eta*power(totalMass,5.0/3.0)*power(piFl,8.0/3.0));  
  tau3 	= pi/(8.0*eta*power(totalMass,2.0/3.0)*power(piFl,5.0/3.0));

  return tau0, tau3


#######################################################
## tauToMasses
#######################################################
def convertTauToMasses( flower, tau0, tau3 ):
  """
  Convert the given masses to the appropriate values of tau0 and tau3.
  Algorithm as above, just reversed
  """
  MTSUN=4.925e-6

  mTotal = 40.0*tau3/( 256.0*pi*pi*tau0*flower)
  eta = 1.0/( 8.0*pow(pi, 2./3.) *pow(flower, 5./3.)*pow(mTotal, 2./3.)*tau3  )

  m1=mTotal*(0.5-sqrt(0.25-eta))
  m2=mTotal-m1

  return m1/MTSUN, m2/MTSUN
  

#######################################################
## calculateFlower
#######################################################
def calculateFlower( m1, m2, tau0 ):
  """
  Calculates the flower involed in the converting
  between tau's and m's.
  Algorithm stolen from 'XLXLALEThincaParameterForInjection'
  from CoincInspiralEllipsoid.cxs
  """
  MTSUN=4.925e-6
  mTotal = m1+m2
  eta=m1*m2/(mTotal*mTotal)
    
  mTotal = mTotal * LAL_MTSUN_SI

  fLower = 5.0 / (256.0 * eta * pow(mTotal, 5.0/3.0) * tau0 )
  fLower = pow(fLower, 3.0/8.0) / pi
  return fLower

##############################################################################
## class grbCoireTable
##############################################################################
class grbTable:
  """
  Table to hold a list of GRB's with all additional informations,
  like time, position, redshift, duration, available segments etc etc.
  """
  class row(object):
    __slots__ = ["identifier", "idRef","gpsTime","rightAscensionDeg","declinationDeg",\
                 "state","type",\
                 "duration", "durationErr","durationRef",\
                 "redshift","redshiftErr","redshiftRef","comments",\
                 "segs","onSource", "onIFOs","siteNumber","onNumber","pmFlag"]
    
    def __init__(self, id="", idRef=-1, gps=-1.0, ra_deg=-1.0, dec_deg=-1.0, state="", type=""):
      """
      Initializing a new GRB event
      """

      # set the identifier
      self.identifier = id
      self.idRef      = idRef

      # set the timeand position
      self.gpsTime    = gps
      self.rightAscensionDeg = ra_deg
      self.declinationDeg    = dec_deg

      # set states and types
      self.state = state # names: POS, REJ, GRB, SGR
      self.type  = type  # short, long

      # set unknown values to each else
      self.duration    = -1.0
      self.durationErr = -1.0      
      self.durationRef = -1
      self.redshift    = -1.0
      self.redshiftErr = -1.0
      self.redshiftRef = -1
      self.comments    = ""

      # segment informations
      self.segs      = None
      self.onSource  = None
      self.onIFOs    = None
      self.onNumber  = 0
      self.siteNumber= 0

    #----------------------------------------------
    # setDuration
    def setDuration(self, dur, durErr, durRef, type = None):

      self.duration    = dur
      self.durationErr = durErr
      self.durationRef = durRef

      if type:
        self.type=type
      else:
        if dur<=2.0:
          self.type="short"
        else:
          self.type="long"
          
    #----------------------------------------------
    # setRedshift
    def setRedshift(self, red, redErr, redRef):
 
      self.redshift    = red
      self.redshiftErr = redErr
      self.redshiftRef = redRef
      
    #----------------------------------------------
    # setOnsource
    def setOnsource(self, left, right):
      """
      Sets the on-source segment
      """
      
      timeOffset = self.gpsTime
      self.onSource = segments.segment( timeOffset-left, timeOffset+right)
      
    #----------------------------------------------
    # setSegments
    def setSegments(self, segs):
      """
      Set the segments and calculates the on-ssource detectors/sites etc
      """
      listIFO=['H1','H2','L1','G1','V1']
      
      onIFOs = ''
      pmFlag=''
      onNumber=0
      
      for ifo in listIFO:
       
        # check if segment overlaps with on-source
        if self.onSource in segs[ifo]:
          onIFOs+=ifo
          onNumber+=1
         
      # determine class
      pmFlag = ''
      if ('H1' in onIFOs) and ('H2' in onIFOs):
        siteNumber=onNumber-1
        pmFlag = '+'
      else:
        siteNumber=onNumber
        if 'H2' in onIFOs:
          pmFlag = '-'

      self.segs      = segs
      self.onIFOs    = onIFOs
      self.onNumber  = onNumber
      self.siteNumber= siteNumber

    #----------------------------------------------
    # setSegments
    def writeExtTrigTable(self, source_file ):
      """
      Creates a XML file containing a ExtTriggersTable with all
      required data.
      """
      
      # prepare a new XML document
      xmldoc = ligolw.Document()
      xmldoc.appendChild(ligolw.LIGO_LW())
      tbl = lsctables.New(lsctables.ExtTriggersTable)
      xmldoc.childNodes[-1].appendChild(tbl)
      
      # get an empty entry and fill it
      row = tbl.RowType()      
      row.process_id    = ''
      row.det_alts      = ''
      row.det_band          = ''
      row.det_fluence       = ''
      row.det_fluence_int   = ''
      row.det_name      = ''
      row.det_peak      = ''
      row.det_peak_int  = ''
      row.det_snr       = ''
      row.email_time    = 0
      row.event_dec     = self.declinationDeg
      row.event_dec_err = 0.0
      row.event_epoch   = ''
      row.event_err_type= ''
      row.event_ra      = self.rightAscensionDeg
      row.event_ra_err  = 0.0
      row.start_time    = int(self.gpsTime)
      row.start_time_ns = int((self.gpsTime - int(self.gpsTime))*1.0e+9)
      row.event_type    = self.type
      row.event_z       = self.redshift
      row.event_z_err   = self.redshiftErr
      row.notice_comments= ''
      row.notice_id     = ''
      row.notice_sequence= ''
      row.notice_time   = 0
      row.notice_type   = ''
      row.notice_url    = ''
      row.obs_fov_dec   = 0.0
      row.obs_fov_dec_width = 0.0
      row.obs_fov_ra    = 0.0
      row.obs_fov_ra_width= 0.0
      row.obs_loc_ele   = 0.0
      row.obs_loc_lat   = 0.0
      row.obs_loc_long  = 0.0
      row.ligo_fave_lho = 0.0
      row.ligo_fave_llo = 0.0
      row.ligo_delay    = 0.0
      row.event_number_gcn = self.idRef
      row.event_number_grb = self.identifier
      row.event_status     = 0
      
      # store the entry in the table and write it to file
      tbl.append(row)
      utils.write_filename(xmldoc, source_file )
      
      
  # ----------------------------------
  # __init__
  def __init__(self, filename = None):
    """
    Initializing the class 'grbTable'  
    """
    self.listIFO=['H1','H2','L1','G1','V1']
    self.rows = []

    self.onsourceLeft = 60
    self.onsourceRight = 120

    if filename:
      self.read(filename)
    
    

  # ----------------------------------
  # write
  def write( self, filename ):
    """
    Write the current table saved to a file
    @param filename: filename which is used
    """

    # open file
    f = open(filename, 'w')
    
    for row in self.rows:
      #line = "%8s %6d %5s %6s %14.3f %8.3f %8.3f  %8.3f %8.3f %6d %7.3f %7.3f %6d %s\n" % \
      #       (row.identifier, row.idRef, row.state, row.type, row.gpsTime,
      #        row.rightAscensionDeg, row.declinationDeg,     
      #        row.duration, row.durationErr, row.durationRef,
      #        row.redshift, row.redshiftErr, row.redshiftRef, row.comments)
      line = "%s,%d,%s,%s,%.0f,%.3f,%.3f,%.3f,%.3f,%d,%.3f,%.3f,%d, %s\n" % \
             (row.identifier, row.idRef, row.state, row.type, row.gpsTime,
              row.rightAscensionDeg, row.declinationDeg,     
              row.duration, row.durationErr, row.durationRef,
              row.redshift, row.redshiftErr, row.redshiftRef, row.comments)
      f.write(line)

    f.close()


  # ----------------------------------
  # read
  def read(self, filename):
    """
    Reads a table from a file
    @param filename: filename which is used
    """
    if False:
      self.rows=[]
      
      # open file
      f = open(filename)
      for lines in f:
        words=lines.split()
        
        # read the specs from the file
        id     = words[0]
        idRef  = int(words[1])
        state  = words[2]
        type   = words[3]
        gps    = float(words[4])
        ra_deg = float(words[5])
        dec_deg= float(words[6])
        dur    = float(words[7])
        durErr = float(words[8])
        durRef = int(words[9])
        z      = float(words[10])
        zErr   = float(words[11])
        zRef   = int(words[12])

        # sanity checks
        if dur>0 and durRef<=0:
          print "WARNING for GRB %s: duration given (%8.3f) but no reference." % (id, dur)
        if z>0 and zRef<=0:
          print "WARNING for GRB %s: redshift given (%8.3f) but no reference." % (id, z)
        
        # define a GRB event
        row = self.row(id, idRef, gps, ra_deg, dec_deg)
        row.setDuration( dur, durErr, durRef, type)
        row.setRedshift( z, zErr, zRef )
        
        # and add it to the list
        self.addRow( row )    

      f.close()

    # read file from xml
    doc = utils.load_filename(filename)    
    extTrigTable = table.get_table(doc, lsctables.ExtTriggersTable.tableName)

    # put in the data
    for entry in extTrigTable:
      gps =  entry.start_time + float(entry.start_time_ns)*1.0e-9
      row = self.row( entry.event_number_grb, entry.event_number_gcn, gps, entry.event_ra, entry.event_dec)
      row.setDuration( -1, 0,0, "unknown")
      row.setRedshift( entry.event_z, entry.event_z_err, 0 )
    
      # and add it to the list
      self.addRow( row ) 
    
    
  # ----------------------------------
  # setSegments
  def setSegments( self, verb=False, dir='SegData' ):
    """
    Read in the segment files and set the segments, on-source detectors etc
    """    
    # loop over each GRB
    for grb in self.rows:
      
      segs={}
      if verb:
        print "Reading segments for GRB%s" % (grb.identifier)

      # set on-source segment
      grb.setOnsource( self.onsourceLeft, self.onsourceRight ) 

      # loop over all IFO's
      onIFOs=''
      onNumber=0
      for ifo in self.listIFO:

        # try to open segment filename
        filename=dir+'/seg'+grb.identifier+'-'+ifo+'.txt'
        try:
          f=open( filename )
        except IOError:
          # either file does not exist, or it is empty. In either way: no segments!
          segs[ifo]=[]
          continue

        # read segment list and close file
        segs[ifo] = segmentsUtils.fromsegwizard(f)
        f.close()

      # set the segment lists
      grb.setSegments( segs )
      if verb:
        print "onIFOs: %s  class: %d-site/%d-detector" % \
              (grb.onIFOs,  grb.siteNumber, grb.onNumber)

    
  # ----------------------------------
  # addRow
  def addRow( self, row ):
    """
    Add a row to the table
    @param row: row to be added
    """

    self.rows.append( row )


  # ----------------------------------
  # find
  def find( self, id ):
    """
    Returns row with the appropriate id
    """
    # loop over all stored GRB's
    for row in self.rows:
      #print "'%s' - '%s'" % (row.identifier, id)
      if row.identifier==id:
        return row

    print "WARNING: GRB "+id+" not found. Returning empty object."
    return None
  



  
##############################################################################
## class grbCoireTable
##############################################################################
class grbCoireTable:
  """
  Table to hold inspiral triggers identified as found/missed injections
  from running with COIRE.
  These are ONE-FILE inputs, created with 'ligolw_add'
  """ 
  
  def __init__(self, foundFile=None, missedFile=None, statistic = None):
    """
    This function reads coincident triggers associated with injections
    (after the COIRE step), assigns a single statistic value to each
    coincidence 
    @param foundFile: name of the file containing the FOUND injections
    @param foundFile: name of the file containing the FOUND injections    
    @param statistic: a structure containing the statistics to use
    """
    self.listIFO=['G1','H1','H2','L1','T1','V1']

    # set the chisq-values
    self.chiThreshold=10
    self.chiBins=16
    self.chiDelta=0.20

    # check the statistic specified
    if not statistic:
      print >>sys.stderr, "No statistic specified, using basic SNR statistic."
      statistic= CoincInspiralUtils.coincStatistic( "snr")
    self.statistic = statistic

    self.found  = []
    self.sim    = []
    self.missed = []

    # vector to mark certain found or missed triggers
    self.markedF = []
    self.markedM = []    
      
    # read found injections if specified
    if foundFile:
      self.addFound( foundFile )

    # read missed injections if specified
    if missedFile:
      self.addMissed( missedFile )

  ## --------------
  ## --addData--
  ## --------------
  def addData( self, foundName, missedName ):
    """
    Adds data to this class.

    @param foundName: name of the COIRE file containing found injections
    @param missedName: name of the file containing missed injections
    """

    # add found injections
    self.addFound( foundName )

    # add missed injections
    self.addMissed( missedName )
    
  ## --------------
  ## --addMissed--
  ## --------------
  def addMissed( self, missedName ):
    """
    Reading the data from a file 'missedName' and adding the data to the
    internal structures.
    @param missedName: name of the file containing missed injections
    """
    # read the file
    missed = self.readMissed( missedName)

    # add the data
    if missed:
      for miss in missed:
        self.missed.append( miss )

        
  ## --------------
  ## --addFound--
  ## --------------
  def addFound( self, foundName ):
    """
    Reading the data from a file 'foundName' and adding the data to the
    internal structures.

    @param foundName: name of the COIRE file containing found injections
    """

    # read the file
    found, sims = self.readFound( foundName)

    # add the data
    if found:
      for coinc in found:
        self.found.append( coinc )

    if sims:
      for sim in sims:
        self.sim.append( sim )

  ## --------------
  ## --readFound--
  ## --------------      
  def readFound( self, foundName ):
    """
    Reads the coincidences from a COIRE file
    @param foundName: name of the COIRE file containing found injections
    """
    # read a xml structure
    doc = utils.load_filename(foundName)
    
    # extract the sim inspiral table
    try: 
      sims = table.get_table(doc, lsctables.SimInspiralTable.tableName)
    except:
      sims = None

    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except:
      snglInspiralTable = None
    
    if ( snglInspiralTable and not sims ) or\
           ( not snglInspiralTable and sims ):
      print >> sys.stderr, "The file "+foundName+\
            " must contain injections AND found triggers!"
      sys.exit(1)

    found = None
    if snglInspiralTable:

      #print snglInspiralTable
      #print len(snglInspiralTable)
      #sys.exit(0)

      # create the coincidences and add the sim-table
      found = CoincInspiralUtils.coincInspiralTable(
        snglInspiralTable, self.statistic )

      # check the equality of sim-tables and trigger-tables
      if len( sims ) != len( found ):
        print >> sys.stderr, "The sim_inspiral table and the coinc_inspiral"\
              " table must have same sizes!"
        sys.exit(1)


      # add the sim tables appropriate
      for index in range( len(sims) ):
        found[index].sim = sims[index]

    return found, sims

  ## --------------
  ## --readMissed--
  ## --------------
  def readMissed( self, missedName ):
    """
    Reads the sim-inspiral tables for the missed injections
    @param missedName: name of the file containing missed injections
    """
    
    # read a xml structure
    doc = utils.load_filename(missedName)
    
    # extract the sim inspiral table
    try: 
      missed = table.get_table(doc, lsctables.SimInspiralTable.tableName)
    except:
      missed = None

    return missed

  ## --------------  
  ## --getSimList--
  ## --------------  
  def getSimList(self, key, index=None, choice="FOUND"):
    """
    Returns a list of injected values, corresponding to the entry with key KEY,
    that are stored in a sim_table in the correct object. 
    Only the indicated entries are returned (INDEX).
    If no INDEX is specified, all values from the whole
    structure/dictionary are returned.
    
    @param key:     name of the field, like snr, mass1 or also totalMass
    @param index:   if specified, returns only values according to this index
    @param choice:  to take the sim triggers either from the FOUND set
                     or the MISSED set
    """

    if choice=="FOUND":
      set=self.sim
    elif choice=="MISSED":
      set=self.missed
    else:
      print >> sys.stderr, "Wrong choice! There are no triggers "+choice+", either FOUND or MISSED"
      sys.exit(1)
      
    if not index:
      index=range(0,len(set))

    list=[]
    if key=="totalMass":
      m1=[getattr(set[i], "mass1") for i in index]
      m2=[getattr(set[i], "mass2") for i in index]      
      list=[m1[i]+m2[i] for i in range(0, len(index))]
      
    elif key=="h_end_time":
      t1=[getattr(set[i], "h_end_time") for i in index]
      t2=[getattr(set[i], "h_end_time_ns") for i in index]
      minTime = min( t1 )
      list=[t1[i]-minTime+t2[i]/1.0e+9 for i in range(0, len(index))]
      
    elif key=="l_end_time":
      t1=[getattr(set[i], "l_end_time") for i in index]
      t2=[getattr(set[i], "l_end_time_ns") for i in index]
      minTime = min( t1 )
      list=[t1[i]-minTime+t2[i]/1.0e+9 for i in range(0, len(index))]
      
    else:
      list=[getattr(set[i], key) for i in index]      

    return list

  
  ## ------------------
  ## ---getTriggerList
  ## ------------------
  def getTriggerList(self, key, ifo, index=None):
    """
    Get a list of values of the found triggers. If 'key' contain the phrase 'end_time',
    then always the end_time related to the specified IFO is used!
    @param key:   Describes the values to return (e.g. 'mass1','totalMass','h_end_time')
    @param ifo:   IFO from which the values should be taken
    @param index: Only the indexed entries are considered. If index=None then all available data will be used. 
    """

    #####################################
    def getList( key, ifo, index ):
      """
      Returns a list of the values wanted
      @param key:   the keyword (like 'mass1','end_time')
      @param ifo:   the detector from which the value should be taken
      @param index: the index of the list-items to be used
      """
      list=[]
      for i in index:
        
        #for coinc in self.found:
        coinc=self.found[i]
        if hasattr(coinc, ifo):
          list.append( getattr(getattr(coinc,ifo),key) )

      return list

    #####################################
    
    if not index:
      index=range(0,len(self.found))

    list=[]
    if key=="totalMass":
      m1=getList( "mass1", ifo, index )
      m2=getList( "mass2", ifo, index )      
      list=[m1[i]+m2[i] for i in range(0, len(m1))]
      
    elif "end_time" in key:
      t1= getList( "end_time" , ifo, index )
      t2= getList( "end_time_ns" , ifo, index )
      minTime = min( t1 )
      list=[t1[i]-minTime+t2[i]/1.0e+9 for i in range(0, len(t1))]  
      
    elif key=="chicut":
      chi = getList( "chisq", ifo, index )
      snr = getList( "snr", ifo, index )      
      threshold=[self.chiThreshold*(self.chiBins+snr[i]*snr[i]*self.chiDelta) for i in range(0, len(chi))]
      list=[chi[i]/threshold[i] for i in range(0, len(chi))]
      
    elif key=="effective_snr":
      for i in index:
        coinc=self.found[i]
        list.append( getattr(coinc,ifo).get_effective_snr() )
        
    else:
      list=getList( key, ifo, index)
    
    return list
  
  ## ----------------
  ## markTriggers
  ## ----------------
  def markTriggers( self, set='missed', parameter='eff_dist_h', order='ascending', number = None):
    """
    Function to mark special triggers, either found/missed triggers and with
    special ordering (ascending/descending). If number is None, all triggers are marked.
    """

    # get the appriopriate data set
    if set=="found":
      valueSet=self.getSimList( parameter, None, "FOUND" )
      list=self.found
    elif set=="missed":
      valueSet=self.getSimList( parameter, None, "MISSED" )
      list=self.missed
    else:
      return

    # sort this list; just get the indices (ascending)
    indexList = argsort( valueSet )

    # reverse order if required
    if order=="descending":
      indexDummy=[indexList[i-1] for i in arange(len(indexList),0,-1)]
      indexList=indexDummy

    # check the 'number' that is set
    if not number:
      number=len(indexList)

    # now populate the indices
    index=[]
    for i in arange(0, number):          
      index.append( indexList[i] )

    # store the index list
    if set=="found":
      self.markedF=index
    elif set=="missed":
      self.markedM=index
    


  ## ------------
  ## sortMasses
  ## ------------
  def sortMasses( self ):
    """
    Function to sort the masses of the found injections, so to ensure that
    mass1 corresponds to the higher mass and mass2 to the lower one
    """

    for i in range(len(self.sim)):

      # retrieve the sim and the coinc entry
      sim=self.sim[i]
      coinc=self.found[i]

      # switch around masses for the simulation, if needed
      if sim.mass2>sim.mass1:
        temp=sim.mass1
        sim.mass1=sim.mass2
        sim.mass2=temp

      # check for any ifo 
      for ifo in self.listIFO:
        if hasattr(coinc, ifo):
          m1= getattr(getattr(coinc,ifo),"mass1")
          m2= getattr(getattr(coinc,ifo),"mass2")

          # switch around the masses, if needed
          if m2>m1:
            setattr(getattr(coinc,ifo),"mass1", m2)
            setattr(getattr(coinc,ifo),"mass2", m1)


            
    
  ## ------------
  ## plotFound
  ## ------------
  def plotFound( self, xvalue='mass1', xifo='H1', xlog='linear',\
                  yvalue='mass2', yifo='H1', ylog='linear',
                  filename=None, flagDiff=None):
    """
    Function to create a plot showing recovered or injected triggers.
    If 'flagDiff=True' the difference between the x and the y values are plotted.
    @param xvalue:  the x-value to be plotted
    @param xifo:    detector of which the data should be plotted. Could also be 'sim'!
    @param xlog:   plotting x-axis linear of logarithmic (linear/log)
    @param yvalue: the y-value to be plotted
    @param yifo:    detector of which the data should be plotted. Could also be 'sim'!
    @param ylog:   plotting y-axis linear of logarithmic (linear/log)
    @param filename: filename to store the plot under
    @param flagDiff: If set to None, the y value is plotted on the y-axis,
                      if set to 'absolute' the absolute difference is plotted (y-x) 
                      and if set to 'relative' the relative difference is plotted ( (y-x)/x )
    """

    ## plot the found injections
    clf()
    if xifo=='sim':
      px=self.getSimList( xvalue )
    else:
      px=self.getTriggerList( xvalue, xifo)
      
    if yifo=='sim':    
      py=self.getSimList( yvalue )
    else:
      py=self.getTriggerList( yvalue, yifo)


    if flagDiff:
      if flagDiff=="absolute":
        pz=[py[i]-px[i] for i in range(0,len(py))]
      if flagDiff=="relative":
        pz=[ (py[i]-px[i])/px[i] for i in range(0,len(py))]
    else:
      pz=py

    # plot the data
    plot( px, pz, 'bo')

    ## insert marked triggers
    if len(self.markedF)>0:
      if xifo=='sim':
        xMarked= self.getSimList( xvalue, self.markedF )
      else:
        xMarked=self.getTriggerList( xvalue, xifo, self.marked )
      if yifo=='sim':   
        yMarked=self.getSimList( yvalue, self.markedF )
      else:
        yMarked=self.getTriggerList( yvalue, yifo, self.marked )

      
      if flagDiff:
        zMarked=[yMarked[i]-xMarked[i] for i in range(0,len(py))]
      else:
        zMarked=yMarked
        
      hold(True)
      plot(xMarked, zMarked, 'mo',markersize=15,\
           markerfacecolor='None', markeredgewidth=2, markeredgecolor='m')
      hold(False)
      
    xlabel( xvalue+' '+xifo, size='x-large')
    if flagDiff:
      ylabel( 'Diff'+flagDiff+': '+yvalue+'_'+yifo+' - '+xvalue+'_'+xifo, size='x-large')
    else:
      ylabel( yvalue+' '+yifo, size='x-large')
    title('Recovered triggers',size='x-large')
    grid(True)
    hold(True)
  
    # scale the axis
    ax=gca()
    ax.set_xscale(xlog)
    ax.set_yscale(ylog)
    
    if filename:      
      savefig(filename+'.png')
    
  ## ----------------------------------
  ## plotScatter
  ## ----------------------------------

  def plotScatter( self, xvalue='mass1', xifo='H1',\
                   yvalue='mass2', yifo='H1', \
                   zvalue='mass2', zifo='H1',\
                   filename=None ):
    """
    Function to create an arbitrary scatter plot
    @param xvalue:   parameters for the x-axis
    @param xifo:     IFO for data on x-axis
    @param yvalue:   parameters for the y-axis
    @param yifo:     IFO for data on y-axis
    @param zvalue:   parameters for the z-axis
    @param zifo:     IFO for data on z-axis    
    @param filename:        name of file 
    """

    # gather the data
    if xifo=='sim':
      px=self.getSimList( xvalue )
    else:
      px=self.getTriggerList( xvalue, xifo)
      
    if yifo=='sim':    
      py=self.getSimList( yvalue )
    else:
      py=self.getTriggerList( yvalue, yifo)

    if zifo=='sim':    
      pz=self.getSimList( zvalue )
    else:
      pz=self.getTriggerList( zvalue, zifo)


    # create plot
    figure(1)
    clf()
    scatter( px, py, 40.0, c=pz, faceted=False )
    colorbar()
    grid(True)
    axis('tight')
    xlabel( xvalue+' '+xifo, size='x-large')
    ylabel( yvalue+' '+yifo, size='x-large')
    title('Scatter plot '+zvalue,size='x-large')
    
    if filename:
      savefig(filename+'.png')
      
  
  ## ----------------------------------
  ## plotScatterCoin
  ## ----------------------------------
  def plotScatterCoin( self, xvalue='mass1', xifo='H1',\
                       yvalue='mass2', yifo='H1', \
                       zvalue='mass2',\
                       filename=None ):
    """
    Function to create an arbitrary scatter plot
    @param xvalue:   parameters for the x-axis
    @param xifo:     IFO for data on x-axis
    @param yvalue:   parameters for the y-axis
    @param yifo:     IFO for data on y-axis
    @param zvalue:   parameters for the z-axis (combined property)
    @param filename:        name of file 
    """

    # gather the data
    if xifo=='sim':
      px=self.getSimList( xvalue )
    else:
      px=self.getTriggerList( xvalue, xifo)

    if yifo=='sim':    
      py=self.getSimList( yvalue )
    else:
      py=self.getTriggerList( yvalue, yifo)

    pz=[]
    if zvalue=='stat':
      for coinc in self.found:
        pz.append( coinc.stat )
        
    elif zvalue=="ethincaHH":
      for coinc in self.found:
        if hasattr( coinc, 'H1') and hasattr( coinc, 'H2'):
          ethinca=tools.XLALCalculateEThincaParameter( coinc.H1, coinc.H2 )
        pz.append( ethinca )
    elif zvalue=="ethincaHL":
      for coinc in self.found:
        if hasattr( coinc, 'H1') and hasattr( coinc, 'L1'):
          ethinca=tools.XLALCalculateEThincaParameter( coinc.H1, coinc.L1 )
        pz.append( ethinca )

    elif zvalue=="edistH1":
      for coinc in self.found:
        if hasattr( coinc, 'H1'):
          edist=tools.XLALEThincaParameterForInjection( coinc.sim, coinc.H1 )
        pz.append( edist ) 
    else:
      pz=None

    # create plot
    clf()
    scatter( px, py, 40.0, c=pz, faceted=False )
    colorbar()
      
    grid(True)
    xlabel( xvalue+' '+xifo, size='x-large')
    ylabel( yvalue+' '+yifo, size='x-large')
    title('Scatter plot '+zvalue,size='x-large')
    
    if filename:
      savefig(filename+'.png')
      
  
  ## -------------------
  ## plotMissedFound
  ## -------------------
  def plotMissedFound( self, xvalue='mass1', xlog='linear',
                       yvalue='mass2', ylog='linear',
                       filename=None, plotLegend=True ):
    """
    Function to create a plot of two parameters from the trigger/sim set,
    with a marker used to distinguish between missed and found.
    The axis can be set to linear or log individually.
    @param xvalue: the x-value to be plotted
    @param xlog:   plotting x-axis linear of logarithmic (linear/log)
    @param yvalue: the y-value to be plotted
    @param ylog:   plotting y-axis linear of logarithmic (linear/log)
    @param filename: filename to store the plot under
    """
   
    ## plot the found injections
    clf()
    xValueFound=self.getSimList( xvalue )    
    yValueFound=self.getSimList( yvalue )
    plot( xValueFound, yValueFound, 'bo')
    xlabel( xvalue, size='x-large')
    ylabel( yvalue, size='x-large')
    title('Missed and found injections',size='x-large')
    grid(True)
    hold(True)
    leg=['found']

    ## insert 'marked found' injections
    if len(self.markedF)>0:
      xValueMarked=self.getSimList( xvalue, self.markedF )
      yValueMarked=self.getSimList( yvalue, self.markedF )      
      plot(xValueMarked, yValueMarked, 'ko',markersize=15,
           markerfacecolor='None', markeredgewidth=2, markeredgecolor='k')
      hold(True)
      leg.append('marked found')

    ## insert 'marked missed' injections
    if len(self.markedM)>0:
      xValueMarked=self.getSimList( xvalue, self.markedM ,"MISSED")
      yValueMarked=self.getSimList( yvalue, self.markedM ,"MISSED" )      
      plot(xValueMarked, yValueMarked, 'go',markersize=15,
           markerfacecolor='None', markeredgewidth=2, markeredgecolor='g')
      hold(True)
      leg.append('marked missed')

    ## plot the missed injections
    if len(self.missed)>0:
      xValueMissed=self.getSimList( xvalue, choice="MISSED" )
      yValueMissed=self.getSimList( yvalue, choice="MISSED" )      
      plot(xValueMissed, yValueMissed, 'rx')
      leg.append('missed')

    # setting the legend
    if plotLegend:
      legend( leg )

    # scale the axis
    ax=gca()
    ax.set_xscale(xlog)
    ax.set_yscale(ylog)     
    if filename:      
      savefig(filename+'.png')

      
  ## -----------------
  ## plotEfficiency
  ## -----------------
  def plotEfficiency( self, value='mass1', nbins=30, range=None, filename=None ):
    """
    Plots the recovery efficiency as a function of the given parameter
    """

    def histng( xdata, xedges, lum_weight=None ):

      """
      histogram xdata with edges specified xedges and yedges.  
      Can rescale the entries by lum_weight
      @param xdata:  array of data for parameter x
      @param xedges: bin boundaries for parameter x
      @param lum_weight: rescaling factor for the histogram
      """
      ng_x = zeros(len(xedges),'d')
      xstep = xedges[1] - xedges[0]

      for i in arange(len(xdata)):
        l = int((xdata[i] - xedges[0])/xstep)

        if not lum_weight:
          lum_array = 1
        else: 
          lum_array = lum_weight[i]
          
        if (l>=0 and l<len(xedges)):
          ng_x[l] += lum_array
 
      return ng_x

    ############################################
    
    valuesFound  = self.getSimList( value, None, choice="FOUND" )    
    valuesMissed = self.getSimList( value, None, choice="MISSED" )
    valuesAll    = self.getSimList( value, None, choice="FOUND" )
    for item in valuesMissed:
      valuesAll.append( item )

    if not range:      
      a=min(valuesAll)
      b=max(valuesAll)
    else:
      a=range[0]
      b=range[1]

    edges=arange(a,b,(b-a)/nbins)      

    # create the data
    all_inj = histng( valuesAll, edges)
    found_inj = histng( valuesFound, edges)
    eff = found_inj / (all_inj + 1e-5)
    
    # create the plot
    clf()
    plot(edges,eff,'b',linewidth=2)

    # include the errors in this plot
    mc = sqrt(found_inj * (all_inj - found_inj) / (all_inj**3 + 1e-5))
    print len(edges), len(eff), len(mc)
    fit= errorbar(edges,eff,mc,fmt=None,ecolor='r',linewidth=2)
    
    xlabel( value, size='x-large')
    ylabel( "efficiency", size='x-large')
    title('Efficiency over '+value,size='x-large')
    grid(True)
    if filename:
      savefig(filename+'.png')
  
   
  ## ------------------
  ## printMissed
  ## ------------------
  def printMissed( self, filename=None, htmlFlag=False):
    """
    Prints out a list of missed injections, unsorted...( by the specified parameter)
    @param filename: file with the results
    @param htmlFlag: if true creates a html table as output
    """

    # open file if required
    file=None
    if filename:
      file=open(filename, 'w')

    # print out
    for i in range(0, len(self.missed)):
      sim=self.missed[i]
      endtime=sim.geocent_end_time+sim.geocent_end_time_ns/1e+9;

      if htmlFlag==False:
        s= "%2d.  distance: %5.1f Mpc |effDistL/H: %5.1f %5.1f Mpc "\
           "|masses: %5.2f  %5.2f |mc: %5.2f |eta: %5.3f |pol: %5.3f  incl: %6.3f "\
           "|time: %6.3f" \
            % (i+1, sim.distance, sim.eff_dist_l, sim.eff_dist_h, sim.mass1, \
               sim.mass2, sim.mchirp, sim.eta, sim.polarization, sim.inclination, endtime)
      else:
        s= "<tr><td> %4d </td><td> %6.2f </td><td> %6.2f </td><td> %6.2f </td><td> "\
           "%6.2f </td><td> %6.2f </td><td> %5.3f </td><td> %6.2f </td><td> "\
           "%6.2f </td><td> </td></tr>" \
            % (i+1, sim.distance, sim.eff_dist_h, sim.mass1, sim.mass2, \
               sim.mchirp, sim.eta, sim.inclination, endtime)

      if file:
        file.write(s+"\n")
      else:
        print s
      
  ## ---------------
  ## printFound
  ## ---------------
  def printFound( self, ifo='H1', filename=None ):
    """
    Function to compile a list of found injections
    """
    c=0

    # open output file if required
    file=None
    if filename:
      file=open(filename, 'w')

    # loop over all found coincidences
    for coinc in self.found:
      
      sim=coinc.sim
      if hasattr( coinc, ifo ):
        trigger=getattr( coinc,ifo )

        c+=1
        endtime=trigger.end_time+trigger.end_time_ns/1e+9;
        s= "%2d. %s |effDist:  %6.2f |masses: %5.2f "\
           "%5.2f sim: %5.2f  %5.2f  |mc: %5.2f sim: %f |eta: %5.3f "\
           "|incl: %5.2f |snr: %6.2f |chisq: %7.2f |effsnr: %6.3f |time: %6.3f" \
             % (c, ifo, trigger.eff_distance, trigger.mass1, \
                trigger.mass2, sim.mass1, sim.mass2, trigger.mchirp, sim.mchirp,\
                trigger.eta, sim.inclination, trigger.snr, trigger.chisq, \
                trigger.get_effective_snr(), endtime)
        if file:
          file.write(s+"\n")
        else:
          print s
      
