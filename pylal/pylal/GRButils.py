import sys
import glob

from pylab import *
from glue.ligolw import utils
from glue.ligolw import table
from glue.ligolw import lsctables
lsctables.SnglInspiralTable.validcolumns.update({
  "ellips_match": "real_8", "wrong": "int_4"
})

import numpy
from pylal import readMeta
from pylal import CoincInspiralUtils
from pylal import SimInspiralUtils
from pylal import SnglInspiralUtils
from pylal import tools
from pylal import support
import PlotUtils

class MySnglInspiral(lsctables.SnglInspiral):
  pass
  __slots__ = lsctables.SnglInspiralTable.validcolumns.keys()

lsctables.SnglInspiralTable.RowType = MySnglInspiral


##-------------------------------------
##-------------------------------------
import string
class NumberReader:
 def __init__(self, file):
  self.fp = file

 def readNumberLine(self):
  line = ''
  while 1:
   line = self.fp.readline()
   if not line:
    break
   if line.strip() == "":
    continue # skip empty lines
   if line.strip().startswith('#'):
    continue # skip comment lines
   return map(float, line.split())

 def readNumbers(self, n=5):
  res = []

  for i in range(n):
   word = ''

   # skip whitespace
   char = self.fp.read(1)
   while char in string.whitespace:
    char = self.fp.read(1)

   # read next word (number)
   while char not in string.whitespace:
    word += char
    char = self.fp.read(1)

   res.append(float(word))

  return res

##############################################################################
# Define a Hist Function
##############################################################################
def histng(xdata,xedges,lum_weight=None):

  """
  histogram xdata with edges specified xedges and yedges.  
  Can rescale the entries by lum_weight
  @param xdata:  array of data for parameter x
  @param xedges: bin boundaries for parameter x
  @param lum_weight: rescaling factor for the histogram
  """
  ng_x = zeros(len(xedges),'d')
  xstep = xedges[1] - xedges[0]
  
  for i in range(len(xdata)):
    l = int((xdata[i] - xedges[0])/xstep)

    if not lum_weight:
      lum_array = 1
    else: 
      lum_array = lum_weight[i]

    if (l>=0 and l<len(xedges)):
      ng_x[l] += lum_array
 
  return ng_x



##############################################################################
## calculate overlap
##############################################################################
def calculateOverlap( trigger1, trigger2, delay=None):
  """
  Calculate the overlap of two sngl_inspiral triggers, given the light-travel time 'delay'.
  Using the equation given in entry
  http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=iags5tdfindchirp&action=view&page=10
  If the standard time-delay is used (Hanford - Livingston) then
  trigger1 must be L and trigger2 must be H!!
  example: calculateOverlap( coinc.L1, coinc.H1, deltaTime )
  @params trigger1: First trigger
  @params trigger2: Second trigger
  @params delay: Time delay between the triggers (in seconds, time2-time1)
  """
  # create the first gamma matrix
  Gamma1=array( [[trigger1.Gamma0, trigger1.Gamma1, trigger1.Gamma2],\
                 [trigger1.Gamma1, trigger1.Gamma3, trigger1.Gamma4],\
                 [trigger1.Gamma2, trigger1.Gamma4, trigger1.Gamma5]])

  # create the second gamma matrix
  Gamma2=array( [[trigger2.Gamma0, trigger2.Gamma1, trigger2.Gamma2],\
                 [trigger2.Gamma1, trigger2.Gamma3, trigger2.Gamma4],\
                 [trigger2.Gamma2, trigger2.Gamma4, trigger2.Gamma5]])

  # create the inverse
  g1=inverse(Gamma1)
  g2=inverse(Gamma2)

  # create the distance vector
  t0=trigger1.tau0-trigger2.tau0
  t3=trigger1.tau3-trigger2.tau3
  tc=trigger1.end_time-trigger2.end_time+(trigger1.end_time_ns-trigger2.end_time_ns)/1.0e+9

  if not delay:
    point=array([tc, t0, t3])
    pointv=reshape( point, (3,1))
  else:
    point=array([tc+delay, t0, t3])
    pointv=reshape( point, (3,1))


  maxVal=0
  # loop over lambda (in 0.01 steps...)
  for lam in arange(0.0, 1.0, 0.01):

    # calculate the 'core' expression
    core=lam*g1+(1-lam)*g2
    icore=inverse(core)
   
    if not delay:
      
      minVal=10000
      for d in arange( -0.01, 0.01, 0.001):
	point=array([tc+d, t0, t3])
        pointv=reshape( point, (3,1))
        overlap= lam*(1-lam)*matrixmultiply( point,matrixmultiply( icore, pointv)  )[0]
        if overlap<minVal:
          minVal=overlap
      overlap=minVal
    else:
      # and the overlap, finally
      overlap= lam*(1-lam)*matrixmultiply( point,matrixmultiply( icore, pointv)  )[0]    

    # find the maximum
    if overlap>maxVal:
      maxVal=overlap

  return maxVal
  
##############################################################################
## calculate range of a binary system 
##############################################################################
def calculateOverlap( mass1=1.4, mass2=1.4, flow=40.0, deltaT=0.000244, deltaF=0.003906, nPoints=1048576, snr=8, psd=None ):
  """
  Calculate the range of a binary system given the PSD and the two m,asses of the system.
  Algorithm taken from the code 'tmpltbank.c'
  @params mass1: Mass of the first component
  @params mass2: Mass of the second component
  @params flower: lower
  @params: psd: 
  """

  # REAL4 compute_candle_distance(REAL4 candleM1, REAL4 candleM2,
  #64 patrick 1.1     REAL4 snr, REAL8 chanDeltaT, INT4 nPoints, 
  #65                 REAL8FrequencySeries *spec, UINT4 cut)
  #66             {
  #67               UINT4 k;
  #68               REAL8 sigmaSqSum = 0;
  #69               REAL8 distance = 0;

  # defining some constants:
  negativeSevenOverThree = -7.0/3.0
  MRSUN=1.4766e+3
  MTSUN=4.925e-6
  PC=3.0857e+16

  distNorm = 2.0 * MRSUN / (1.0e6 * PC )
  
  # what else do I need?: deltaT and deltaF and nPoints
  # cut=int(flower/deltaF)
  # numPoint=2*length(psd), might not...
  
  totalMass =  mass1 + mass2
  mu = mass1 * mass2 / totalMass
  a = sqrt( (5.0 * mu) / 96.0 ) * pow( totalMass/( pi*pi ), 1.0/3.0 ) * pow( MTSUN / deltaT, -1.0/6.0 )
  sigmaSq = 4.0 * ( deltaT / nPoints ) * distNorm * distNorm * a * a
  fmax = 1.0 / (6.0 * sqrt(6.0) * pi * totalMass * MTSUN)
  f = 0
  
  sigmaSqSum = 0
  binLow=int(flow/deltaF)
  binHigh=int(fmax/deltaF)
  print binLow, binHigh
  for k in range(binLow, binHigh):
    f=k*deltaF
    sigmaSqSum+=pow( float(k) /  nPoints, negativeSevenOverThree ) / psd[k]

  sigmaSq *= sigmaSqSum;
  distance = sqrt( sigmaSq ) / snr;

  return distance
 
 #82               for ( k = cut, f = spec->deltaF * cut; 
 #83                   k < spec->data->length && f < fmax; 
 #84                   ++k, f = spec->deltaF * k )
 #85 patrick 1.1   {
 #86                 sigmaSqSum += 
 #87                   pow( (REAL8) k / (REAL8) nPoints, negativeSevenOverThree ) 
 #88                   / spec->data->data[k];
 #89               }


##############################################################################
## class coincGRBInspiralTable
##############################################################################
class coincGRBInspiralTable(CoincInspiralUtils.coincInspiralTable):
  """
  Table to hold one set of coincidence triggers, as read from one THINCA or COIRE file
  """

  def __init__(self,  inspTriggers = None, stat = None):
    """
    Initializing this table
    @param inspTriggers: a metaDataTable containing inspiral triggers 
                         from which to construct coincidences
    @param stat:         an instance of coincStatistic
    """
    self.listIFO=['G1','H1','H2','L1','T1','V1']
    self.slideTimes= {'H1': 0.0, 'H2': 10.0, 'L1':5.0}

    # initialize the underlying class (the original coincInspiralTable)
    CoincInspiralUtils.coincInspiralTable.__init__(self,inspTriggers, stat)

    #self.reslide(821428403,600)

    # number for accounting the number of the injection
    self.injNumber=None

    # flag to indicate if this coincidence was 'found'
    self.found=False

    ## calculate the sllipse-match between each trigger and the sim
    #self.calcEllipsMatch()

  ## ---reslide---##
  def reslide(self, gpsStartTime, gpsLength):
    """
    reslides triggers.
    NEED TO SLIDE THEM ON THE RING!!!
    """
    for coinc in self.rows:
      idNumber=coinc.event_id
      numSlide=(( idNumber % 1000000000) / 100000)
      if numSlide>5000:
        numSlide=5000-numSlide

      # loop through all possible IFO's (H1, H2 and L1 only!!!)
      for ifo,slideT in self.slideTimes.iteritems():
		
      	if hasattr(coinc, ifo):
          slideTime=numSlide*slideT
          trigger=getattr( coinc, ifo )
          newTime=trigger.end_time+slideTime

          # check the time
          if newTime<gpsStartTime:
            newTime+=gpsLength
          if newTime>gpsStartTime+gpsLength:
            newTime-=gpsLength

          # restore the time
          trigger.end_time=newTime
          setattr( coinc, ifo, trigger )   

      #if hasattr(coinc, 'H1') and hasattr(coinc, 'L1'):
         #print idNumber, numSlide, coinc.H1.end_time- coinc.L1.end_time

   ## ---getTriggerListIntern--- ##
  def getTriggerListIntern(self, key, ifo):
    """
    Get a list of values from the stored coincidences. The values related to 'key'
    are returned for an IFO 'ifo'. In the normal case only one value should be returned
    (in the case of clustered data), but is also SHOULD works else.
    If index=None then all available data will be used. 
    """
    values=[]
    for coinc in self.rows:
      if hasattr(coinc, ifo):
        values.append( getattr(getattr(coinc,ifo),key) )

    return values
        
    
  ## ---setInjNumber ---##
  def setInjNumber(self, number):
    """
    Sets the actual number of injection, i.e. the number associated with the filename
    and the random seed
    @param number: number of injection
    """    
    self.injNumber=number

  ## ---getInjNumber ---##
  def getInjNumber(self):
    """
    Returns the actual number of injection, i.e. the number associated with the filename
    and the random seed
    """
    return self.injNumber

  ## ---setSimTable ---##
  def setSimTable(self, sim_inspiral):
    """
    Sets the sim_inspiral table associated with these coincidences. Many coincidences allowed with one sim table
    @param sim_inspiral: sim_inspiral table to set
    """
    self.sim_table=sim_inspiral
    
  ## ---getSimTable ---##
  def getSimTable(self):
    """
    Returns the current sim_inspiral table
    """
    return self.sim_table
  
  ## --- setFound ---##
  def setFound( self ):
    """
    Sets FOUND flag
    """
    self.found=True
  
  ## --- setMissed ---##
  def setMissed( self ):
    """
    Sets MISSED flag
    """
    self.found=False

   ## --- isFound ---##
  def isFound( self ):
    """
    Returns FOUND flag
    """
    return self.found

  ## --- calcEllipsMatch---##
  def calcEllipsMatch( self ):
    """
    Calculating the ellipsoid match between each coincidence and the sim table
    """
    
    # loop over each interferometer
    for ifo in self.listIFO:
  
      # loop over each coincident trigger
      for coinc in self.rows:
        if hasattr(coinc, ifo):

          # extract the trigger and the value
          trigger=getattr(coinc,ifo)
          sim=self.sim_table

          # check the input values for safety
          if type(trigger)!=MySnglInspiral:
            print >>sys.stderr, "Wrong type for 'trigger' in coincGRBInspiralTable.calcEllipsMatch: "+str(type(trigger))
            sys.exit(1)
          if type(sim)!=lsctables.SimInspiral:
            print >>sys.stderr, "Wrong type for 'sim' in coincGRBInspiralTable.calcEllipsMatch: "+str(type(trigger))
            sys.exit(1)
                   
          # call the function from LAL to calculate the overlap
          par=tools.XLALEThincaParameterForInjection(sim, trigger)

          # store this information
          trigger.ellips_match=par

  

##############################################################################
# function to read in a list of files and extract the simInspiral tables
##############################################################################
def readFiles(fileGlob,statistic=None):
  """
  read in the Sngl and SimInspiralTables from a list of files
  if Sngls are found, construct coincs, add injections (if any)
  also return Sims (if any)
  @param fileGlob: glob of input files
  @param statistic: statistic to use in creating coincs
  """
  #if fileGlob is empty return empty structures...
  if not fileGlob:
    if opts.verbose:
      print "Warning: No glob specified, returning empty structures..."
    return None, coincGRBInspiralTable() 

  fList = glob.glob(fileGlob)
  if not fList:
    print >>sys.stderr, "The glob for " + fileGlob + " returned no files"
    sys.exit(1)
  sims = None
  coincs = None
  for thisFile in fList:

    # read a xml structure
    doc = utils.load_filename(thisFile)

    # extract the sim inspiral table
    try: 
      simInspiralTable = \
          table.get_table(doc, lsctables.SimInspiralTable.tableName)
      sims = simInspiralTable
    except: simInspiralTable = None


    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: snglInspiralTable = None
    if snglInspiralTable:

      # create the coincidences and add the sim-table
      coincInspiralTable = \
        coincGRBInspiralTable(snglInspiralTable,statistic)
      if simInspiralTable:
        for i in range(len(coincInspiralTable)):
          coincInspiralTable[i].add_sim(simInspiralTable)

    else:
      # no inspiral triggers found, initialize structure anyway and
      # add missed sim
      coincInspiralTable= coincGRBInspiralTable()
      if simInspiralTable: 
        coincInspiralTable.add_missed_sims(simInspiralTable)

    # add the coincidences to 'coincs'
    if coincs:
      coincs.extend(coincInspiralTable)
    else:
      coincs = coincInspiralTable
          
  return sims,coincs, snglInspiralTable

      
##############################################################################
## class singleInjections
##############################################################################
class singleInjections:
  """
  Table to hold single injection triggers from FOUND and MISSED xml files
  """
  
  ## ----------------------------------
  ## __init__
  ## ----------------------------------
  def __init__(self, dir, userTag=None, injNumbers = None, multiTriggers = False):
    """
    Initialise the structure by reading the H1/H2/L1 files from the specified directory.
    A user-tag can be specified, so the following files are being read:
      ifo+'-SIRE-'+userTag+str(inj)+'_FOUND.xml
    The value 'inj' loops through the injections I want to read.
    If 'multiCoire' is set to true then multiple triggers are allowed for each file.
    @params dir:     directory that holds the files
    @params userTag: optional usertag
    @params injNumber:  specifies the injections to be read
    @params multiCoire: allows multiple triggers if true
    """

    # create dictionaries
    self.dataDict={}
    self.dataDict["H1"]={}
    self.dataDict["H2"]={}
    self.dataDict["L1"]={}        
    
    self.ifoList=['H1','H2','L1']

    # store the input values
    self.injNumbers=injNumbers

    # loop over all three IFO's
    for ifo in self.ifoList:

      # create basic file name
      baseName=dir+ifo+'-SIRE-'
      if userTag:
        baseName=baseName+userTag

      foundInj=[]
      missedInj=[]

      # loop over the injection range
      injList=range(injNumbers[0], injNumbers[1])
      for inj in injList:

        help={}
        help["sngl"]=[]
        help["sim"]=[]
                
        # create filename and read data
        filename=[baseName+str(inj)+'_FOUND.xml']
        sim1=SimInspiralUtils.ReadSimInspiralFromFiles(filename)
        sngl=SnglInspiralUtils.ReadSnglInspiralFromFiles(filename)

        text='Injection '+str(inj)+' for detector '+ifo+' ... ' # just a comment

        # store sngl table                
        if sngl:
          
          # check if more than one trigger is being read
          if not multiTriggers and len(sngl)>1:
            print "ERROR: More than ONE inspiral trigger in file "+filename[0]
            sys.exit(1)

          help["sngl"]=sngl
          text=text+'found'    # just a comment
          foundInj.append( inj )
          
        else:
          text=text+'missed'  # just a comment 
          missedInj.append( inj )
        print text   # just a comment, print it

        # create filename and read data in case the injection is missed
        filename=[baseName+str(inj)+'_MISSED.xml']
        sim2=SimInspiralUtils.ReadSimInspiralFromFiles(filename)

        # store sim table (whereever this table is from)
        if not sim2:
          help["sim"]=sim1[0]
        else:
          help["sim"]=sim2[0]

        
        # put the data into the dictionary
        self.dataDict[ifo][inj]=help

      self.dataDict[ifo]["found"] = foundInj
      self.dataDict[ifo]["missed"]= missedInj
      self.dataDict[ifo]["wrong"]= []
      
    # find the wrong injections, do this by users' request
    #self.checkWrongInjections()
      
  ## ----------------------------------
  ## checkWrongInjections
  ## ----------------------------------
  def checkWrongInjections(self, matchCut=50):
    """
    Checking the stored triggers if there is a 'wrong' found one.
    @params matchCut: the cut value
    """
    
    # loop over all three IFO's
    for ifo in self.ifoList:

      # re-evaluate found and missed
      foundInj=[]
      missedInj=[]
      wrongInj=[]
      prevFound=self.dataDict[ifo]["found"]
      
      for inj in range(self.injNumbers[0], self.injNumbers[1]):

        sngl_table=self.dataDict[ifo][inj]["sngl"]
        sim_table=self.dataDict[ifo][inj]["sim"]

        foundFlag=False
        for trigger in sngl_table:

          # calculate the match between sim and trigger...
          match=tools.XLALEThincaParameterForInjection( sim_table, trigger )
          trigger.ellips_match=match
          if match>matchCut:
            trigger.wrong=1
          else:
            foundFlag=True
            trigger.wrong=0          

        if foundFlag:
          foundInj.append(inj)
          print "Injection "+str(inj)+ "for detector "+ifo +" is found"
        else:
          missedInj.append(inj)
          print "Injection "+str(inj)+ "for detector "+ifo +" is missed"
          if inj in prevFound:
            wrongInj.append(inj)

      # store new 'found' and 'missed'    
      self.dataDict[ifo]["found"] = foundInj
      self.dataDict[ifo]["missed"]= missedInj
      self.dataDict[ifo]["wrong"]= wrongInj      

  
  ## ----------------------------------
  ## plotMissedFound
  ## ----------------------------------
  def plotMissedFound( self, xvalue='mchirp', linlog='lin',userDist='eff',filename=[]):
    """
    Function to create a plot of found/missed injections for all IFO's
    @params xvalue: value to be plotted on the x-axis. Every valid sim_imspiral entry
                   AND totalMass (m1+m2) can be plotted
    @params linlog: specifies the kind of plot (linear/logarithmic)
    @params userDist: specifies the distance to use. either the effective distance (eff) or the real distance( real)
    @params filename: filename to save plot in. If not specified no plot will be saved
    """

    # loop over all three IFO's
    for ifo in self.ifoList:

      # check userDist
      if userDist=="real":
        useDist="distance"
        nameDist="Real distance"
      else:
        if ifo=="L1":        
          useDist="eff_dist_l"
          nameDist="Effective Distance L"
        else:
          useDist="eff_dist_h"
          nameDist="Effective Distance H"

      # generate list of found and missed injections
      listFound = self.dataDict[ifo]["found"]
      listMissed= self.dataDict[ifo]["missed"]
      listWrong=  self.dataDict[ifo]["wrong"]

      # create non-intersecting list for the missed injections
      dummy=[]
      for e in listMissed:
        if e not in listWrong:
          dummy.append(e)
      listMissed=dummy
      
      # now get the data to plot, first the y-values (distance)
      ylist0=[ getattr(self.dataDict[ifo][i]["sim"],useDist) for i in listMissed ]
      ylist1=[ getattr(self.dataDict[ifo][i]["sim"],useDist) for i in listWrong ]      
      ylist2=[ getattr(self.dataDict[ifo][i]["sim"],useDist) for i in listFound ]

      # now the x-values
      if xvalue=="totalMass": 
        xlist0=self.getValList( self.dataDict[ifo], "sim", xvalue, listMissed)
        xlist1=self.getValList( self.dataDict[ifo], "sim", xvalue, listWrong)
        xlist2=self.getValList( self.dataDict[ifo], "sim", xvalue, listFound)                   
      elif xvalue=="end_time":
        if ifo=="L1":
	  xval="l"
        else:
          xval="h"
        xlist0=self.getValList( self.dataDict[ifo], "sim", xval, listMissed)
        xlist1=self.getValList( self.dataDict[ifo], "sim", xval, listWrong)
        xlist2=self.getValList( self.dataDict[ifo], "sim", xval, listFound)
      else:
        xlist0=[ self.dataDict[ifo][i]["sim"].__getattribute__(xvalue) for i in listMissed ]
        xlist1=[ self.dataDict[ifo][i]["sim"].__getattribute__(xvalue) for i in listWrong ]
        xlist2=[ self.dataDict[ifo][i]["sim"].__getattribute__(xvalue) for i in listFound ]
                          
      # create the plot
      clf()
      lw=5
      if linlog=='lin':
        semilogy( xlist0, ylist0,'r^',linewidth=lw)
        hold(True)
        semilogy( xlist2, ylist2, 'bo')        
        semilogy( xlist1, ylist1, 'gD')        
      else:
        loglog( xlist0, ylist0,'r^',linewidth=lw)
        hold(True)
        loglog( xlist2, ylist2, 'bo')
        loglog( xlist1, ylist1, 'gD' )     

      hold(False)
      xlabel(xvalue, size='x-large')
      ylabel(nameDist,size='x-large')
      grid(True)
      legend(('missed','found','wrong'))
      title('Missed and found injections '+ifo, size='x-large' )
      if filename:
        savefig(filename+'_'+xvalue+'_'+ifo+'.png')

        
  ## ----------------------------------
  ## plotAccuracy
  ## ----------------------------------
  def plotAccuracy( self, nbins, itemX, itemY=None, mode="rel", cluster=False, filename=''):
    """
    Function to create accuracy plots of a property between injected and recovered
    @params nbins: number of bins for the histogram plot
    @params itemX: item of which to plot the accuracy
    @params itemY: optional parameter to use as a y dimension
    @params mode:  if 'rel' plots the relative accuracy, if 'abs' plots the absolute difference
    @poarams cluster: if True the triggers are being clustered accordong to the SNR
    @filename:     name of file 
    """

    # loop over all three IFO's
    for ifo in self.ifoList:

      ifoPrefix=ifo[0].lower()
      
      # generate list of found and missed injections
      listFound = self.dataDict[ifo]["found"]
      listMissed= self.dataDict[ifo]["missed"]
      listWrong=  self.dataDict[ifo]["wrong"]

      # create non-intersecting list for the missed injections
      dummy=[]
      for e in listMissed:
        if e not in listWrong:
          dummy.append(e)
      listMissed=dummy

      simsFound =[ self.dataDict[ifo][i]["sim"]  for i in listFound ]
      snglsFound=[ self.dataDict[ifo][i]["sngl"] for i in listFound ]
      arrFound=arange(0,len(simsFound))

      simsWrong =[ self.dataDict[ifo][i]["sim"]  for i in listWrong ]
      snglsWrong=[ self.dataDict[ifo][i]["sngl"] for i in listWrong ]
      arrWrong=arange(0,len(simsWrong))

      dataX, dataY=self.getAccuracyValues(arrFound, simsFound, snglsFound, ifo, itemX, itemY, cluster, mode)
      wrongX, wrongY=self.getAccuracyValues(arrWrong, simsWrong, snglsWrong, ifo, itemX, itemY, cluster, mode)      
    
      # do the plot here
      if len(dataX)>0:

        clf()
        if wrongX:
          minVal=min( [min(dataX), min(wrongX)])
          maxVal=min( [max(dataX), max(wrongX)])
          width=(maxVal-minVal)/nbins
          xVector=arange( minVal, maxVal, width)
          y0, x0,dummy=hist(dataX, xVector)
          y1, x1,dummy=hist( wrongX, xVector)
          bar( x0, y0, color='b',width=0.8*width)
          hold(True)
          x1=x1+(x1[1]-x1[0])*0.2
          bar( x1, y1, color='g',width=0.8*width)
          hold(False)
        else:
          hist(dataX, nbins)
                
        xlabel('Delta '+itemX+' '+ifo)
        ylabel('#')
        title('Accuracy plot '+itemX)
        grid(True)
        if filename:
          savefig(filename+'_'+itemX+'_'+ifo+'.png')
          
        # create a SNR dependence plot
        if len(dataY)>0:
          clf()

          plot( dataX, dataY, 'bo')
          if wrongX:
            hold(True)
            plot( wrongX, wrongY, 'gD')
            hold(False)
          xlabel('Delta '+itemX+' '+ifo)
          ylabel(itemY+' '+ifo)
          title('Accuracy plot '+itemX)
          grid(True)
          if filename:
            savefig(filename+'_'+itemX+'-'+itemY+'_'+ifo+'.png')
          
  ## ----------------------------------
  ## plotAccuracy
  ## ----------------------------------
  def getAccuracyValues( self, arrFound, simsList, snglsList, ifo, itemX, itemY, cluster, mode):
    """
    Return values used to plot for accuracies
    """
    ifoPrefix=ifo[0].lower()
    # prepare lists to hold the actual data
    dataX=[]
    dataY=[]

    for index in arrFound:
      
      # get the sim and sngl table(s)
      sim=simsList[index]
      sngl_table=snglsList[index]
      
      # retrieve the simulated value
      if itemX=="end_time":
        injected=getattr( sim, ifoPrefix+'_end_time')+\
                  getattr( sim, ifoPrefix+'_end_time_ns')/1e+9
      elif itemX=="mtotal":
        injected=sim.mass1+sim.mass2
      else:
        injected=getattr( sim, itemX)

      if cluster:
        triggerList=[self.getCluster( sngl_table )]
      else:
        triggerList=sngl_table

      # loop over all triggers (might be one only...)
      for trigger in triggerList:

        if not trigger:
          continue
          
        # retrieve the recovered and injected values
        if itemX=='end_time':
          recovered=trigger.end_time+trigger.end_time_ns/1e+9;
        elif itemX=="mtotal":
          recovered=trigger.mass1+trigger.mass2
        else:
          recovered=getattr(trigger, itemX)
            
        if mode=="rel":
          # add relative difference
          dataX.append((recovered-injected)/injected)
        else:
          # or absolute difference
          dataX.append(recovered-injected)

        if itemY:
          dataY.append( getattr(trigger, itemY) )
            
    return dataX, dataY
    
  ## ----------------------------------
  ## getValList
  ## ----------------------------------
  def getValList(self, object, table, key, index):
    """
    INTERNAL function returns a list of values (labeled 'key') from a 'table' of an 'object'. Only the
    values which are in a 'index' list are returned.    
    """
    list=[]
    if key=="totalMass":
      m1=[object[i][table].__getattribute__("mass1") for i in index]
      m2=[object[i][table].__getattribute__("mass2") for i in index]
      list=[m1[i]+m2[i] for i in range(0, len(index))]
    elif key=="l" or key=="h":
      list=[object[i][table].__getattribute__(key+"_end_time")+object[i][table].__getattribute__(key+"_end_time_ns")/1.0e-9 for i in index]
    else:
      list=[object[i][table].__getattribute__(key) for i in index]

    return list

  ## ----------------------------------
  ## getCluster
  ## ----------------------------------
  def getCluster(self, sngl_list, index=None):
    """
    Return the cluster-trigger from a list of triggers
    @params sngl_list: list of triggers
    @params index: 
    """
    # create index
    if not index:
      index=arange(1, len( sngl_list))

    maxSNR=0
    maxIndex=None
    # loop over all triggers to look at
    for i in index:

      # check for new maximum
      if sngl_list[i].snr>maxSNR:
        maxSNR=sngl_list[i].snr
        maxIndex=i

    if maxIndex:
      return sngl_list[maxIndex]
    else:
      return None
  


class TestClass(CoincInspiralUtils.coincInspiralTable.row):
  __slots__ = ["event_id", "numifos","stat","G1","H1","H2",\
                 "L1","T1","V1","sim","rsq","bl","inj","found","wrong"]

  
##############################################################################
## class grbFullTable
##############################################################################
class grbFullTable(CoincInspiralUtils.coincInspiralTable):
  """
  Table to hold inspiral triggers obtained by a GRB run of the pipeline.
  """ 

  ## ----------------------------------
  ## __init__
  ## ----------------------------------
  def __init__(self, globFiles=None, statistic = None, num_slides=50, veto_file=None):
    """
    @param glob: a glob containing all triggers (full data and timeslides)
    @param stat: a structure containing the statistics to use
    """

    #class row(CoincInspiralUtils.coincInspiralTable.row):
    #  __slots__ += ["inj","found","wrong"]
    
    #CoincInspiralUtils.coincInspiralTable.__init__(self,None, statistic)
    self.num_slides=num_slides
    self.rows=[]

    # check if glob is specified
    if globFiles:

      # read the file(s)
      sims,coincs,singles=readFiles( globFiles, statistic )
      for row in coincs.rows:
        self.rows.append(row)
      #self.append(coincs)

      # initialize thr mother class to store the singles
      CoincInspiralUtils.coincInspiralTable.__init__(self,singles, statistic)
    else:
      CoincInspiralUtils.coincInspiralTable.__init__(self,None, statistic)

      
      print len(coincs)

  ## ----------------------------------
  ## applyVeto
  ## ----------------------------------
  def applyVeto(self, vetoFile):
    """
    Applying a veto to the data and finding coincs again
    """
    if vetoFile:
      file = open( vetoFile , "r")
      seglist = segmentsUtils.fromsegwizard(file)
      file.close()
      self.inspTriggers = self.inspTriggers.veto(seglist,"end_time")
      self.inspSlides   = self.inspSlides.veto(seglist,"end_time")

    # refind coincident triggers
    getCoincidences(self.statistic)

  ## ----------------------------------
  ## cluster
  ## ----------------------------------
  def cluster( self,clusterTime ):
    """
    Cluster the coincident triggers
    """
    self.coincTriggers = self.coincTriggers.cluster( clusterTime )
    self.coincSlides   = self.coincSlides.cluster( clusterTime )

  ## ----------------------------------
  ## getCoincidences
  ## ----------------------------------
  def getCoincidences(self):
    """
    Find the coincident triggers
    """
    # analyze zero-lag triggers
    if self.inspTriggers:
      self.coincTriggers = readMeta.coincInspiralTable(self.inspTriggers, self.statistic)
      print  self.coincTriggers

    # analyze time-slide triggers
    if self.inspSlides:
      slide_num = range(1 , self.num_slides + 1)
      slide_num.extend(range(-self.num_slides, 0))

      for slide in slide_num:
        this_slide = {}
        this_slide["slide_num"] = slide
        if slide > 0:
          this_slide["sngl_trigs"] = self.inspSlide.getslide(slide)
        else:
          this_slide["sngl_trigs"] = self.inspSlide.getslide(5000 - slide)

      # make coincs
      this_slide["coinc_trigs"] = \
        readMeta.coincInspiralTable(this_slide["sngl_trigs"], self.statistic )

      # add slide to list
      self.coincSlides.append(this_slide)

  ## ----------------------------------
  ## plotChisquare
  ## ----------------------------------
  def plotChisquare( self, ifo, filename):
    """
    Function to create a plot of found/missed injections
    """    
    data=self.getsngls(ifo)
    snr=  [ d.snr for d in data]
    chisq=[ d.chisq for d in data]
    
    p=PlotUtils.LogyPlot( snr, chisq,'bo')
    p.setLabels('SNR','Chisq','x-large')
    p.setTitle('SNR and Chisquare values')
    grid(True)
    p.draw()
  

    if filename:
      p.savepic(filename+'.png')     


##############################################################################
## class grbInjectionTable
##############################################################################
class grbInjectionTable:
  """
  Table to hold coincident inspiral triggers and the corresponding injection triggers.
  This can be either from THINCA files or from COIRE files.
  """ 
  
  def __init__(self, globName=None, statistic = None, injNumbers = None, globSim = None, multiCoire = False):
    """
    This function reads coincident triggers associated with injections
    (after the COIRE step), assigns a single statistic value to each
    coincidence and checks for wrong injections.
    @param glob: a glob for the FOUND/MISSED files
    @param statistic: a structure containing the statistics to use
    @param injNumbers: 2-component vector to specify the injections to read.
    """
    self.listIFO=['G1','H1','H2','L1','T1','V1']

    # set the chisq-values
    self.chiThreshold=10
    self.chiBins=16
    self.chiDelta=0.20
    
    if not injNumbers:
      print >>sys.stderr, "No injection numbers specified. Must specify range of injections to be read."
      sys.exit(1)

    if not statistic:
      print >>sys.stderr, "No statistic specified, using basic SNR statistic."
      statistic= CoincInspiralUtils.coincStatistic( "snr")

    # dummy storage
    self.stat=statistic
    self.injNumbers=injNumbers

    # main list holding the data
    self.list=[]

    # lists containing the 'indices' of the found, missed and wrong-found coincidences
    # 'found' and 'missed' are exclusive sets, 'wrong' is a subset of 'found' 
    self.found=[]
    self.wrong=[]
    self.missed=[]
      
    # read the data
    if globSim:
      self.readThinca(globName, globSim)
    else:
      self.readCoire(globName, multiCoire)

  ## --------------
  ## --readThinca--
  ## --------------      
  def readThinca( self, globName, globSim):
    """
    Reads the coincidences from THINCA files and the corresponding injection (HL) files.
    @param globName: glob for the THINCA files
    @param simFiles: glob for the HL files
    """

    # loop over the injections to read
    for inj in range(self.injNumbers[0], self.injNumbers[1]):

      # read the THINCA-file
      filename=globName+"_"+str(inj)+"-*"
      print filename
      dummy,coinc, dummy = readFiles(filename, self.stat)    
        
      # read the HL-files
      filename=globSim+"."+str(inj)+".xml"
      print filename
      sim,dummy,dummy = readFiles(filename, self.stat)

      if len(sim)>1:
        print >> sys.stderr, "sim_inspiral table with more than one sim read from file "\
              +filename+". Not allowed."
        sys.exit(1)

      # append new object to the list
      coinc.setSimTable( sim[0] )
      coinc.setInjNumber( inj )
      coinc.calcEllipsMatch()  # calculate the ellipse-match
      self.list.append( coinc )

  ## --------------
  ## --readCoire--
  ## --------------
  def readCoire( self, globName, multiCoire):
    """
    Reads the coincident triggers from the COIRE files, identify the missed and found injections,
    checks for wrong injections TO BE DONE LATER!!!)
    @params globName: glob for the missed and found injections
    @params multiCoire: if set to false, at maximum ONE coincidence expected per COIRE file
    """

    if globName:

      self.rows=[]
      # loop over the injections
      for inj in range(self.injNumbers[0],self.injNumbers[1]):

        print "Reading COIRE "+str(inj)+"..."
        
        # read the found-files
        filename=globName+"-"+str(inj)+"_FOUND.xml"
        simF,coincF,singlesF=readFiles(filename, self.stat)
        
        # read the missed-files
        filename=globName+"-"+str(inj)+"_MISSED.xml"
        simM,coincM,singlesM=readFiles(filename, self.stat)
          
        # coincidences found
        coinc=None
        if coincF:
          coinc=coincF
          sim=simF
          coinc.setFound()
        elif coincM:
          coinc=coincM
          sim=simM
          coinc.setMissed()

        if coinc:
          # append sim table and append whole structure to 'itself'
          if len(coinc.rows)>1 and not multiCoire:
            print "WARNING: "+str(len(coinc.rows))+\
                  " coincidences (rows) found for injection "+str(inj)

          # error if more than one sim table found
          if len(sim)>1:
            print "ERROR: "+str(len(sim))+\
                " injections found for injection "+str(inj)
            sys.exit(1)

          coinc.setSimTable( sim[0] )
          coinc.setInjNumber( inj )
          coinc.calcEllipsMatch()
          
          # add the coincInspiralTable to the master list
          self.list.append(coinc)
        else:
          print "No COINC found"
          
    # finally, loop over ALL read triggers, identify the
    # 'found' and 'missed' injections
    for i in range( 0,len(self.list) ):
      if self.list[i].isFound():
        self.found.append(i)
        print "injection #"+str(self.list[i].getInjNumber())+" is found"
      else:
        self.missed.append(i)
        print "injection #"+str(self.list[i].getInjNumber())+" is missed"
      
    # output to screen
    print "Triggers found: "+str(len(self.found))+ " triggers missed: "+\
          str(len(self.missed))

  ## --------------
  ## --getTimeDiff--
  ## --------------
  def getTimeDiff(self, trigger, sim):
    """
    Returns the difference of time between the given trigger and
    the given injection in seconds. Uses the correct IFO in sim.
    Returnvalue = simTime - triggerTime
    @trigger: structure containing one sngl_inspiral trigger
    @sim    : structure containing one sim_inspiral trigger
    """

    if trigger.ifo=="L1":
      simTime  = sim.l_end_time
      simTimeNS= sim.l_end_time_ns
    if trigger.ifo=="H1" or trigger.ifo=="H2" :
      simTime  = sim.h_end_time
      simTimeNS= sim.h_end_time_ns

    diff= (trigger.end_time- simTime)+(trigger.end_time_ns-simTimeNS)/1.0e+9
    return -diff

  ## --------------
  ## --getValList--
  ## --------------
  def getValList(self, object, key, index=None):
    """
    Returns a list of values that are stored in a
    structure/dictionary OBJECT and has a kay name KEY.
    Only the indicated entries are returned (INDEX).
    If no INDEX is specified, all values from the whole
    structure/dictionary are returned.

    @object:  list of structures/dictionary, as a sngl_inspiral
              or sim_inspiral table
    @key:     name of the field, like snr, mass1 or also totalMass
    @index:   if specified, returns only values according to this index
    """
    if not index:
      index=range(0,len(object))

    list=[]
    if key=="totalMass":
      m1=[object[i].__getattribute__("mass1") for i in index]
      m2=[object[i].__getattribute__("mass2") for i in index]
      list=[m1[i]+m2[i] for i in range(0, len(index))] 
    else:
      list=[object[i].__getattribute__(key) for i in index]

    return list

  ## ------------------
  ## ---getTriggerList
  ## ------------------
  def getTriggerList(self, key, ifo, index=None):
    """
    Get a list of values from the recovered triggers. 'key' describes what values to return.
    If 'key; contain the phrase 'end_time', then always the end_time related to the specified IFO is used!!
    If index=None then all available data will be used. 
    """
    from operator import add
    
    if not index:
      index=range(0,len(self.list))

    list=[]
    if key=="totalMass":
      m1=[self.list[i].getTriggerListIntern( "mass1", ifo)  for i in index]
      m2=[self.list[i].getTriggerListIntern( "mass2", ifo)  for i in index]
      m1=reduce(add, m1)
      m2=reduce(add, m2)
      list=[m1[i]+m2[i] for i in range(0, len(m1))]
    elif "end_time" in key:
      t1=[self.list[i].getTriggerListIntern( "end_time", ifo)  for i in index]
      t2=[self.list[i].getTriggerListIntern( "end_time_ns", ifo)  for i in index]
      t1=reduce(add, t1)
      t2=reduce(add, t2)
      list=[t1[i]+t2[i]/1.0e+9 for i in range(0, len(t1))]
    elif key=="chicut":
      chi=[self.list[i].getTriggerListIntern( "chisq", ifo)  for i in index]
      snr=[self.list[i].getTriggerListIntern( "snr", ifo)  for i in index]
      chi=reduce(add, chi)
      snr=reduce(add, snr)
      threshold=[self.chiThreshold*(self.chiBins+snr[i]*snr[i]*self.chiDelta) for i in range(0, len(chi))]
      list=[chi[i]/threshold[i] for i in range(0, len(chi))]
    elif key=="effective_snr":
      for i in index:
        for coinc in self.list[i].rows:
          if hasattr( coinc,ifo ):
            list.append( getattr(coinc,ifo).get_effective_snr() ) 
      #list=[self.list[i].get_effective_snr()  for i in index]
    else:
      list=[self.list[i].getTriggerListIntern( key, ifo)  for i in index]
      list=reduce(add, list)
    
    return list

  ## --------------  
  ## --getSimList--
  ## --------------  
  def getSimList(self, key, index=None):
    """
    Returns a list of values, corresponding to the entry with key KEY,
    that are stored in a sim_table in the correct object. 
    Only the indicated entries are returned (INDEX).
    If no INDEX is specified, all values from the whole
    structure/dictionary are returned.
    
    @key:     name of the field, like snr, mass1 or also totalMass
    @index:   if specified, returns only values according to this index
    """
    if not index:
      index=range(0,len(self.list))

    list=[]
    if key=="totalMass":
      m1=[getattr(self.list[i].getSimTable(), "mass1") for i in index]
      m2=[getattr(self.list[i].getSimTable(), "mass2") for i in index]      
      list=[m1[i]+m2[i] for i in range(0, len(index))]
    elif key=="h_end_time":
      t1=[getattr(self.list[i].getSimTable(), "h_end_time") for i in index]
      t2=[getattr(self.list[i].getSimTable(), "h_end_time_ns") for i in index]
      list=[t1[i]+t2[i]/1.0e+9 for i in range(0, len(index))]
    elif key=="l_end_time":
      t1=[getattr(self.list[i].getSimTable(), "l_end_time") for i in index]
      t2=[getattr(self.list[i].getSimTable(), "l_end_time_ns") for i in index]
      list=[t1[i]+t2[i]/1.0e+9 for i in range(0, len(index))]
    else:
      list=[getattr(self.list[i].getSimTable(), key) for i in index]      

    return list
    
  
  ## -----------------
  ## plotMissedFound
  ## -----------------
  def plotMissedFound( self, xvalue='mchirp', linlog='lin',
                       useDist='distance', filename=None ):
    """
    Function to create a plot of found/missed injections
    """
   
    ## plot the found injections
    clf()
    simDistF=[getattr(self.list[i].getSimTable(), useDist) for i in self.found]
    xValueFound=self.getSimList( xvalue, self.found )    
    if linlog=='lin':
      semilogy(xValueFound, simDistF, 'bo')
    else:
      loglog( xValueFound, simDistF, 'bo')
    xlabel( xvalue, size='x-large')
    ylabel('Distance',size='x-large')
    title('Missed and found injections',size='x-large')
    grid(True)
    hold(True)

    ## insert 'wrong' injections
    if len(self.wrong)>0:
      simDistW=[getattr(self.list[i].getSimTable(), useDist) \
                for i in self.wrong]
      xValueWrong=self.getSimList( xvalue, self.wrong )
      if linlog=='lin':
        semilogy(xValueWrong, simDistW, 'mo',markersize=15,\
                 markerfacecolor=None, markeredgewidth=2, markeredgecolor='m')
      else:
        loglog( xValueWrong, simDistW,'mo',markersize=15,\
                markerfacecolor=None, markeredgewidth=2, markeredgecolor='m')
      hold(True)

    ## plot the missed injections
    if len(self.missed)>0:
      simDistM=[getattr(self.list[i].getSimTable(), useDist) \
                for i in self.missed]
      xValueMissed=self.getSimList( xvalue, self.missed )
      if linlog=='lin':
        semilogy(xValueMissed, simDistM, 'rx')
      else:
        loglog( xValueMissed, simDistM, 'rx')

    plotLine=0
    if plotLine:
      # calculate a line that DISTinguish the found/missed regime
      lineX, lineY = self.computeDistLine( xValueFound, simDistF, \
                                           xValueMissed, simDistM)
      lineX[0]=lineX[0]+0.001
      hold(True)
      semilogy( lineX, lineY, 'g-', linewidth=3)
      hold(False)
    
    if len(self.wrong)>0:
      legend( ('found','marked','missed') )
    else:
      legend( ('found','missed') )

    if xvalue=='h_end_time':
      axis([854378304,854378844, 1, 10000 ]) 
    if filename:      
      savefig(filename+'.png')
      
  ## -------------------
  ## plotMissedFound2
  ## -------------------
  def plotMissedFound2( self, xvalue='mass1', xlog='linear', yvalue='mass2', ylog='linear',
                        filename=None ):
    """
    Function to create a plot of two parameters from the trigger/sim set,
    with a marker used to distinguish between missed and found.
    The axis can be set to linear or log individually.
    """
   
    ## plot the found injections
    clf()
    xValueFound=self.getSimList( xvalue, self.found )    
    yValueFound=self.getSimList( yvalue, self.found )
    plot( xValueFound, yValueFound, 'bo')
    xlabel( xvalue, size='x-large')
    ylabel( yvalue, size='x-large')
    title('Missed and found injections',size='x-large')
    grid(True)
    hold(True)

    ## insert 'wrong' injections
    if len(self.wrong)>0:
      xValueWrong=self.getSimList( xvalue, self.wrong )
      yValueWrong=self.getSimList( yvalue, self.wrong )      
      plot(xValueWrong, yValueWrong, 'mo',markersize=15,markerfacecolor=None, markeredgewidth=2, markeredgecolor='m')
      hold(True)

    ## plot the missed injections
    if len(self.missed)>0:
      xValueMissed=self.getSimList( xvalue, self.missed )
      yValueMissed=self.getSimList( yvalue, self.missed )      
      plot(xValueMissed, yValueMissed, 'rx')
       
    if len(self.wrong)>0:
      legend( ('found','marked','missed') )
    else:
      legend( ('found','missed') )

    #if xvalue=='h_end_time':
    #  axis([853168884,853169064, 1, 10000 ])

    # scale the axis
    ax=gca()
    ax.set_xscale(xlog)
    ax.set_yscale(ylog)     
    if filename:      
      savefig(filename+'.png')

      
  ## ------------
  ## plotFound
  ## ------------
  def plotFound( self, xvalue='mass1', xifo='H1', xlog='linear', yvalue='mass2', yifo='H1', ylog='linear',
                        filename=None ):
    """
    Function to create a plots showing recovered triggers.
    For plots showing sim/triggers see plotSimFound
    """

    ## plot the found injections
    clf()
    px=self.getTriggerList( xvalue, xifo)
    py=self.getTriggerList( yvalue, yifo)
    plot( px, py, 'bo')

    ## insert marked triggers
    if len(self.wrong)>0:
      xMarked=self.getTriggerList( xvalue, xifo, self.wrong )
      yMarked=self.getTriggerList( yvalue, yifo, self.wrong )
      hold(True)
      plot(xMarked, yMarked, 'mo',markersize=15,\
           markerfacecolor=None, markeredgewidth=2, markeredgecolor='m')
      hold(False)
      
    xlabel( xvalue+' '+xifo, size='x-large')
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

  ## ------------
  ## plotFound2
  ## ------------
  def plotFound2( self, xvalue='mass1', xifo='H1', xlog='linear',\
                  yvalue='mass2', yifo='H1', ylog='linear',
                        filename=None ):
    """
    Function to create a plots showing recovered triggers.
    For plots showing sim/triggers see plotSimFound
    """

    ## plot the found injections
    clf()
    if xifo=='sim':
      px=self.getSimList( xvalue, self.found)
    else:
      px=self.getTriggerList( xvalue, xifo)
    if yifo=='sim':    
      py=self.getSimList( yvalue, self.found)
    else:
      py=self.getTriggerList( yvalue, yifo)


    # plot the data
    plot( px, py, 'bo')

    ## insert marked triggers
    if len(self.wrong)>0:
      if xifo=='sim':
        xMarked= self.getSimList( xvalue, self.wrong)
      else:
        xMarked=self.getTriggerList( xvalue, xifo, self.wrong )
      if yifo=='sim':   
        yMarked=self.getSimList( yvalue, self.wrong)
      else:
        yMarked=self.getTriggerList( yvalue, yifo, self.wrong )
        
      hold(True)
      plot(xMarked, yMarked, 'mo',markersize=15,\
           markerfacecolor=None, markeredgewidth=2, markeredgecolor='m')
      hold(False)
      
    xlabel( xvalue+' '+xifo, size='x-large')
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


  ## ---------------
  ## plotSimFound
  ## ---------------
  def plotSimFound( self, xvalue='mass1', xlog='linear', yvalue='mass2', yifo='H1', ylog='linear',
                        filename=None ):
    """
    Function to create a plots showing injected and recovered values for a certain IFO.
    """
  
    ## plot the found injections
    clf()
    px=self.getSimList( xvalue, self.found )   
    #px=self.getTriggerList( xvalue, xifo)
    py=self.getTriggerList( yvalue, yifo)
    diff=[py[i]-px[i] for i in range(0,len(py))]
    plot( px, diff, 'bo')


    ## insert marked triggers
    if len(self.wrong)>0:
      xMarked=self.getSimList( xvalue, self.wrong )
      yMarked=self.getTriggerList( yvalue, yifo, self.wrong )
      print xMarked, yMarked
      diffMarked=[yMarked[i]-xMarked[i] for i in range(0,len(yMarked))]
      hold(True)
      plot(xMarked, diffMarked, 'mo',markersize=15,\
           markerfacecolor=None, markeredgewidth=1, markeredgecolor='m')
      hold(False)
    
    xlabel( xvalue+' inj', size='x-large')
    ylabel( yvalue+' '+yifo+' -' + xvalue+' inj', size='x-large')
    title('Injected and recovered triggers',size='x-large')
    grid(True)
    hold(True)
  
    # scale the axis
    ax=gca()
    ax.set_xscale(xlog)
    ax.set_yscale(ylog)
    axis('tight')
    if filename:      
      savefig(filename+'.png')
      
  ## -----------------
  ## plotEfficiency
  ## -----------------
  def plotEfficiency( self, value='mass1', nbins=30, range=None, filename=None ):
    """
    Plots the recovery efficiency as a function of the given parameter
    """
    valuesFound  = self.getSimList( value, self.found )    
    valuesMissed = self.getSimList( value, self.missed )
    valuesAll    = self.getSimList( value)

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
    fit,errors = errorbar(edges,eff,mc,fmt=None,ecolor='r',linewidth=2)
    
    xlabel( value, size='x-large')
    ylabel( "efficiency", size='x-large')
    title('Efficiency over '+value,size='x-large')
    grid(True)
    if filename:
      savefig(filename+'.png')
  
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
      valueSet=self.getSimList( parameter, self.found )
      list=self.found
    elif set=="missed":
      valueSet=self.getSimList( parameter, self.missed )
      list=self.missed
    else:
      return

    #print valueSet
    #print list

    # sort this list; just get the indices (ascending)
    indexList = numpy.argsort( valueSet )
    #print indexList

    # reverse order if required
    if order=="descending":
      indexDummy=[indexList[i-1] for i in arange(len(indexList),0,-1)]
      indexList=indexDummy
    #print indexList

    # check the 'number' that is set
    if not number:
      number=len(indexList)

    # now create the 'special'-list
    for i in arange(0, number):
      self.wrong.append( list[indexList[i]] )
      #print  list[indexList[i]]

    # test:
    #valueSet=self.getSimList( parameter, self.wrong )
    #print valueSet

    
  ## -----------------
  ## computeDistLine
  ## -----------------    
  def computeDistLine2( self, xFound, yFound, xMissed, yMissed):
    """
    returns coordinates for a line that DISTinguish the found/missed injections
    @xFound:  x value of found injections
    @yFound:  y value of found injections
    @xMissed: x value of missed injections
    @yMissed: y value of missed injections
    """

    # other possible algorithms see: /home/alex/Work/E_ExtTrigger/Temp/testerMiddle

    # retrieve min/max for the x-axis
    xMin=min( [min(xFound), min(xMissed)] )
    xMax=max( [max(xFound), max(xMissed)] )
    n=20 # divide into 20 'bins' for plotting the line
    sigmasq=10
    xDelta=(xMax-xMin)/n

    # create list containing the points for the lines
    xArray=[]
    yArray=[]

    # loop over the number of bins
    for i in range(0,n):
      xPoint=i*xDelta+xMin



      testMin=xPoint-xDelta
      testMax=xPoint+xDelta

      # find the distances relevant to the x-axis area
      vecMissed=[yMissed[j] for j in range(len(yMissed))
                 if xMissed[j]>testMin and xMissed[j]<testMax ]
      vecFound= [yFound[j] for j in range(len(yFound))
                 if xFound[j]>testMin and xFound[j]<testMax ]
                 
      # store this point
      if vecMissed or vecFound:

        if vecMissed:
          minMissed=min( vecMissed)
          dist=minMissed
        if vecFound:
          maxFound =max( vecFound )
          dist=maxFound

        if vecMissed and vecFound:
          dist=(minMissed+maxFound)/2.0
      
        xArray.append(xPoint+xDelta/2.0)
        yArray.append(dist)

    #print yArray
    sigma=200/(xMax-xMin)
    temp=copy.deepcopy(yArray)
    for i in range(0,len(xArray)):
      xpoint=xArray[i]

      yArray[i]=0
      norm=0;
      for j in range(0,len(xArray)):
        delta=xpoint-xArray[j]
        ampli=exp(-delta*delta/sigma)
        yArray[i]+=ampli*temp[j]        
        norm+=ampli
        
      yArray[i]=yArray[i]/norm
    #print yArray

    return xArray, yArray

  ## -------------------
  ## computeDistLine
  ## -------------------
  def computeDistLine( self, xFound, yFound, xMissed, yMissed):
    """
    returns coordinates for a line that DISTinguish the found/missed injections
    @xFound:  x value of found injections
    @yFound:  y value of found injections
    @xMissed: x value of missed injections
    @yMissed: y value of missed injections
    """

    # other possible algorithms see: /home/alex/Work/E_ExtTrigger/Temp/testerMiddle

    # retrieve min/max for the x-axis
    xMin=min( [min(xFound), min(xMissed)] )
    xMax=max( [max(xFound), max(xMissed)] )
    n=20   # divide into 20 'bins' for plotting the line
    #sigmasq=10
    xDelta=(xMax-xMin)/n

    # create list containing the points for the lines
    xArray=[]
    yArray=[]

    # loop over the number of bins
    for i in range(0,n):
      xPoint=i*xDelta+xMin

      # set the test points
      testMin=xPoint-2*xDelta
      testMax=xPoint+2*xDelta
      
      # find the distances relevant to the x-axis area
      rangeMissed=[yMissed[j] for j in range(len(yMissed))
                 if xMissed[j]>testMin and xMissed[j]<testMax ]
      rangeFound= [yFound[j] for j in range(len(yFound))
                 if xFound[j]>testMin and xFound[j]<testMax ]

      # set value to zero
      R=0
      for dist in range(1, 500, 5):
        vecMissed=[ rangeMissed[j]for j in range(len(rangeMissed))
                 if rangeMissed[j]<=dist]
        vecFound=[ rangeFound[j]for j in range(len(rangeFound))
                 if rangeFound[j]<=dist]

        nm=len(vecMissed)
        nf=len(vecFound)
        sum=nm+nf
        if sum>0:
          ratio=float(nf)/float(sum)
        else:
          ratio=1.0

        if ratio<=0.9 and R==0:
          R=dist
        print "bin %d,  dist: %d  numbers: %d/ %d dist: %f  R:%f" \
              % (i, dist, len(vecMissed), len(vecFound), ratio, R)
     
      xArray.append(xPoint+xDelta/2.0)
      yArray.append(R)
      print "bin %d,  dist: %d" % (i, R)
    #sys.exit(0)

    return xArray, yArray
      
    
  ## ------------------
  ## printMissed
  ## ------------------
  def printMissed( self, filename=None ):
    """
    Prints out a list of missed injections, sorted by the distance
    @filename: file with the results
    """

    file=None
    if filename:
      file=open(filename, 'w')

    #index=self.indexMissed
    simDistM=[getattr(self.list[i].getSimTable(), "distance") for i in self.missed]
    index=argsort(simDistM)
    sims=[self.list[i].getSimTable() for i in take(self.missed, index)]
    injNumber=[self.list[i].getInjNumber() for i in take(self.missed, index)]
     
    # print out
    for i in range(0, len(sims)):
      sim=sims[i]
      inj=injNumber[i]
      #number=sim.inj#self.missed[index[i]]
      endtime=sim.geocent_end_time+sim.geocent_end_time_ns/1e+9;
      s= "%2d. ( %3d ) distance: %6.2f Mpc | effDist L/H: %6.2f / %6.2f | masses: %5.2f / %5.2f mc: %5.2f | eta: %5.3f | pol: %5.3f  incl: %6.3f | time: %6.3f" \
            % (i+1,inj,sim.distance,sim.eff_dist_l, sim.eff_dist_h, sim.mass1, sim.mass2, sim.mchirp, sim.eta, sim.polarization, sim.inclination, endtime)

      if file:
        file.write(s+"\n")
      else:
        print s
      
  ## ---------------
  ## printFound
  ## ---------------
  def printFound( self, ifo='H1', index=None, filename=None ):
    """
    Function to compile a list of found injections
    """
    c=0

    file=None
    if filename:
      file=open(filename, 'w')

    if not index:
      index=range(0, len(self.list))

    for i in index:
      inj=self.list[i]
      injNumber=self.list.index(inj)
      sim=inj.getSimTable()
      for coinc in inj:
        if hasattr(coinc, ifo):
          trigger=getattr(coinc,ifo)

          c+=1
          endtime=trigger.end_time+trigger.end_time_ns/1e+9;
          s= "%2d. ( %3d ) %s | effDist:  %6.2f | masses: %5.2f / %5.2f ( %5.2f / %5.2f ) mc: %5.2f ( %f )| eta: %5.3f | snr: %6.2f chisq: %7.2f | effsnr: %6.3f | time: %6.3f" \
             % (c,injNumber, ifo, trigger.eff_distance, trigger.mass1, trigger.mass2, sim.mass1, sim.mass2, trigger.mchirp, sim.mchirp,trigger.eta, trigger.snr, trigger.chisq, trigger.get_effective_snr(), endtime)
          if file:
            file.write(s+"\n")
          else:
            print s
      
    
##############################################################################
## class grbInjectionTable
##############################################################################
class grbAllTable(grbInjectionTable):
  """
  Table to hold inspiral injection triggers obtained by a GRB run of the pipeline,
  as well as the full-data(plgr) and background triggers
  """
  ## ----------------------------------
  ## __init__
  ## ----------------------------------
  def __init__(self, globInj=None, globFull=None, globSlide=None, statistic = None, injNumbers = None):
    """
    This function reads coincident trigger from injections, zero-lag and timeslides.
    It used the given statistic to assign a single statistic value to each coincident trigger. 
    @param globInj: a glob for the FOUND/MISSED files
    @param globFull: a glob for the full (plgr) files
    @param globSlide: a glob for the background files    
    @param statistic: a structure containing the statistics to use
    """
    # initialize the injection part
    grbInjectionTable.__init__(self, globInj, statistic, injNumbers)

    self.listIFO=['G1','H1','H2','L1','T1','V1']
    self.listIFOtags=['g', 'h', 'h', 'l', 't', 'v']
    
    # variable declaration
    self.coincFG=[]
    self.singlesFG=None
    self.coincBG=[]
    self.singlesBG=None
    
    # initialize the other parts
    if globFull:
      self.initFull(globFull)
    if globSlide:
      self.initBackground(globSlide)
      
  ## ----------------------------------
  ## initInjections
  ## ----------------------------------
  def initInjections(self, globInj=None):
    """
    @globInj: glob for the injections
    Reading the injections from a directory
    """
    self.initInj( self, globInj)
    
  ## ----------------------------------
  ## initFull
  ## ----------------------------------
  def initFull( self, globFull=None):
    """
    @globFull: glob for the full data triggers (or playground)
    Reading the full data (playground) from a directory
    """
    # check if glob is specified
    if globFull:

      # read the file(s)
      sims,self.coincFG,self.singlesFG=readFiles( globFull, self.stat )
      print str(len(self.coincFG)) + " foreground triggers read"

     
  ## ----------------------------------
  ## initBackground
  ## ----------------------------------
  def initBackground( self, globBG=None):
    """
    @globBG: glob for the background triggers
    Reading the background triggers
    """
    # check if glob is specified
    if globBG:

      # read the file(s)
      sims,self.coincBG,self.singlesBG=readFiles( globBG, self.stat )
      print str(len(self.coincBG)) + " background triggers read"

  ## ----------------------------------
  ## getSingles 
  ## ----------------------------------
  def getSingles(self, coincList, ifo=None, ifo2=None ):
    """
    return the single-IFO values to create e.g. scatter plots
    @coincList:  list of coincidences
    @ifo:        interferometer; if not specified all ifo's are used
    """
    singleList=[]
    for coinc in coincList:

      if ifo2:
        if hasattr(coinc, ifo) and hasattr(coinc, ifo2):
          singleList.append(getattr(coinc,ifo))
      else:

        if ifo:
          if hasattr(coinc, ifo):
            singleList.append(getattr(coinc,ifo))
        else:
          for ifo in self.listOFO:
            if hasattr(coinc, ifo):
              singleList.append(getattr(coinc,ifo))

    return singleList
      
  ## ----------------------------------
  ## plotChisquare
  ## ----------------------------------
  def plotChisquare( self, ifo, filename):
    """
    Function to create a plot of found/missed injections
    """
    # create dummy plot
    p=PlotUtils.Plot( 0, 0,'w.')
    p.setLabels('SNR','Chisq','x-large')
    p.setTitle('SNR and Chisquare values')
    p.drawLabels()
    hold(True)

    import pickle
    
    # draw the different items
    if len(self.coincBG):
      data=self.getSingles( self.coincBG, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'kx')
      #f=open('bg_snr.dat', 'w')
      #f.write(str(snr))
      #f.close()
      #f=open('bg_chi.dat', 'w')
      #f.write(str(chisq))
      #f.close()
      
    if len(self.coincFG):
      data=self.getSingles( self.coincFG, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'ro')
      #f=open('fg_snr.dat', 'w')
      #f.write(str(snr))
      #f.close()
      #f=open('fg_chi.dat', 'w')
      #f.write(str(chisq))
      #f.close()

    if len(self.rows):
      data=self.getSingles( self.rows, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'bo')
      #f=open('inj_snr.dat', 'w')
      #f.write(str(snr))
      #f.close()
      #f=open('inj_chi.dat', 'w')
      #f.write(str(chisq))
      #f.close()
      
    grid(True)
    hold(False)
    #axis([0,0, 25, 200])
    xlim(0,25)
    ylim(0,200)

    if filename:
      p.savepic(filename+'.png')

  ## ----------------------------------
  ## plotTimeSeries
  ## ----------------------------------
  def plotTimeSeries( self, dataset, ifo, variable, filename=None):
    """
    Function to plot a time series of some variables for the given
    ifo and the given dataset
    @dataset:  full, slides or injections
    @ifo:      the IFO
    @variable: specifies what to plot (e.g. snr)
    @filename: filename to save the picture in (optional)
    """

    data=None
    if dataset=='full':
      data=self.getSingles( self.coincFG, ifo)
    elif dataset=='slides':
      data=self.getSingles( self.coincBG, ifo)
    elif dataset=='injections':
      data=self.getSingles( self.rows, ifo)

    if data:
      time= [ d.end_time+d.end_time_ns/1.0e+9 for d in data]
      yval= [ d.__getattribute__(variable) for d in data]
      p=PlotUtils.Plot( 0, 0,'w.')
      p.setLabels('time',variable,'x-large')
      p.setTitle('not ready yet...')
      p.drawLabels()
      p.draw( time, yval,'ro')

      if filename:
        p.savepic(filename+'_timeSeries.png')  

    
  ## ----------------------------------
  ## plotScatterSNR
  ## ----------------------------------
  def plotScatterSNR( self, ifo, ifo2, filename):
    """
    Function to create a plot of found/missed injections
    """
    # create dummy plot
    p=PlotUtils.Plot( 0, 0,'w.')
    p.setLabels('SNR','Chisq','x-large')
    p.setTitle('SNR and Chisquare values')
    p.drawLabels()
    hold(True)

    import pickle
    
    # draw the different items
    if len(self.coincBG):
      data=self.getSingles( self.coincBG, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'kx')
      
    if len(self.coincFG):
      data=self.getSingles( self.coincFG, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'ro')

    if len(self.rows):
      data=self.getSingles( self.rows, ifo)
      snr=  [ d.snr for d in data]
      chisq=[ d.chisq for d in data]
      p.draw( snr, chisq,'bo')
      
    grid(True)
    hold(False)
    #axis([0,0, 25, 200])
    xlim(0,25)
    ylim(0,200)

    if filename:
      p.savepic(filename+'_snr_snr.png')     
      
  ## ----------------------------------
  ## plotScatter 
  ## ----------------------------------
  def plotScatter( self, xaxis, yaxis, title, filename):
    """
    Function to create an arbitrary scatter plot
    @xaxis:   parameters for the x-axis
    @yaxis:   parameters for the y-axis    
    @title:   title of the plot
    @filename: name of file 
    """
    # create dummy plot
    clf()
    p=PlotUtils.Plot( 0, 0,'w.')
    p.setLabels(xaxis.getLabel(),yaxis.getLabel(),'x-large')
    p.setTitle(title)
    p.drawLabels()
    hold(True)

    #import pickle
    
    # draw the different items
    coincSet=[self.coincBG, self.coincFG, self.rows]
    colorSet=['kx','ro','bo']
    for c in range(0,3):
      if coincSet[c]:
        #print c, len(coincSet[c])
          
        datax=self.getSingles( coincSet[c], xaxis.ifo, yaxis.ifo)
        datay=self.getSingles( coincSet[c], yaxis.ifo, xaxis.ifo)      
        xval= [ d.__getattribute__(xaxis.name) for d in datax]
        yval= [ d.__getattribute__(yaxis.name) for d in datay]
        p.draw( xval, yval, colorSet[c])
        #if yaxis.ifo=='H2' and c==0:
        #  print len(datax)
        #  print xval
        #  print yval
        #  sys.exit(0)

        #if c==1:
        #  print xval, yval
          
    grid(True)
    hold(False)

    # set the axis
    if xaxis.range:
      xlim( xaxis.range[0], xaxis.range[1] )
    if yaxis.range:
      ylim( yaxis.range[0], yaxis.range[1] )

    if filename:
      p.savepic(filename+'_'+xaxis.name+'_'+yaxis.name+'.png')
      
  ## ----------------------------------
  ## plotAccuracies
  ## ----------------------------------
  def plotAccuracies( self, nbins, filename):
    """
    Function to plot graphs showing the accuracy between injected
    and recovered signals for dist, masses, time, ...
    @filename: name of file 
    """
  
    # loop over each item
    itemList=['end_time','mchirp','eta']
    for item in itemList:
      plotAccuracy( nbins, item, 'snr', 0, filename)
      

  ## ----------------------------------
  ## plotAccuracy
  ## ----------------------------------
  def plotAccuracy( self, nbins, item, itemX, rel=0, filename=''):
    """
    Function to plot one graphs showing the accuracy between one
    injected and recovered prpperty
    @filename: name of file 
    """
    # get list of coincs
    coincList=[self.rows[i] for i in self.found]
    
    # loop over each interferometer
    for ifo in self.listIFO:

      vecX=[]
      vecDiff=[]

      # loop over each coincident trigger
      for coinc in coincList:
        if hasattr(coinc, ifo):

          # extract the trigger and the value
          trigger=getattr(coinc,ifo)
          
          # retrieve the recovered and injected values
          if item=='end_time':
            recovered=trigger.end_time+trigger.end_time_ns/1e+9;
            i = self.listIFO.index(ifo)
            itemSim = self.listIFOtags[i]+'_'
            injected=getattr(coinc.sim, itemSim+'end_time')+\
                      getattr(coinc.sim, itemSim+'end_time_ns')/1e+9
          else:
            recovered=getattr(trigger, item)
            injected=getattr(coinc.sim, item)
            
          if rel==1:
            # add relative difference
            vecDiff.append((recovered-injected)/injected)
          else:
            # or absolute difference
            vecDiff.append(recovered-injected)
          vecX.append( getattr(trigger, itemX) )
            
            
      # do the plot here
      if len(vecDiff)>0:

        clf()
        print vecDiff
        hist(vecDiff, nbins)
        xlabel(item+' '+ifo)
        ylabel('#')
        title('Accuracy plot '+item)
        grid(True)
        if filename:
          savefig(filename+'_'+item+'_'+ifo+'.png')
          
        # create a SNR dependence plot
        clf()
        print vecX
        print vecDiff
        sys.exit()
        p=PlotUtils.Plot(vecX, vecDiff,'bo')
        p.setLabels(itemX,'Delta '+item+' '+ifo)
        p.setTitle('Accuracy plot '+item)
        p.draw()
        if filename:
          p.savepic(filename+'_'+item+'-'+itemX+'_'+ifo+'.png')
          
    
    ## ----------------------------------
  ## 
  ## ----------------------------------

    ## ----------------------------------
  ## 
  ## ----------------------------------
