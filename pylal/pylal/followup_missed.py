# $Id$
#
# Copyright (C) 2006  Alexander Dietz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#__Id__ = "$Id$"
#__author__ = "Darren Woods and Stephen Fairhurst <sfairhurs@gravity.phys.uwm.edu>"
#__version__ = "$Revision$"[11:-2]
#__date__ = "$Date$"[7:-2]
#__name__ = "followup_missed.py"
#__title__ = "Followup missed injections"

import os, sys, exceptions, copy
from math import sqrt, pi

import matplotlib
matplotlib.use('Agg')
from pylab import rcParams, fill, figtext, figure, plot, axes, axis, xlabel, ylabel, title, close, grid, legend
try:
  set
except NameError:
  from sets import Set as set

from pylal import SnglInspiralUtils
from pylal import SimInspiralUtils
from pylal import CoincInspiralUtils
from pylal import InspiralUtils
from glue import lal
from glue import markup
from glue import segments
from glue import segmentsUtils
from glue.markup import oneliner as extra
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
import numpy


##########################################################
class FollowupMissed:
  """
  This defines a class for followup missed injections.
  The definition as a class interfers minimal with the original code.

  Usage:

  # first initialize the class
  followup = followup_missed.FollowupMissed( cache, miscache, opts)

  # later, one can do a followup of a missed injection 'inj', which returns the filename
  # of the html file containing the tables and pictures created
  followuphtml = followup.followup( inj, ifo )
  """


  # -----------------------------------------------------
  def __init__(self, cache, coireVetoMissedCache, opts, \
               deltaTime=20,injectionWindow = 50):
    """
    Initialize this class and sets up all the cache files.
    @param cache: the cache of all files
    @param coireVetoMissedCache: the cache that contains the COIRE-files
                                  of the missed injections
    @param opts: The 'opts' from the main code
    @param deltaTime: The time window that is used in the time-series plots (in seconds)
    """
    rcParams.update({'text.usetex': False})

    self.colors = {'H1':'r','H2':'b','L1':'g','V1':'m','G1':'c'}
    self.stageLabels = ['TMPLTBANK', 'INSPIRAL_FIRST','THINCA_FIRST',\
                        'TRIGBANK', 'INSPIRAL_SECOND', 'THINCA_SECOND']
    self.orderLabels = copy.deepcopy(self.stageLabels)
    self.orderLabels.extend( [ 'THINCA_SECOND_CAT_1','THINCA_SECOND_CAT_2', \
                               'THINCA_SECOND_CAT_3','THINCA_SECOND_CAT_4'] )

    # setting all the parameters
    self.cache                = cache
    self.coireVetoMissedCache = coireVetoMissedCache
    self.opts                 = opts
    self.deltaTime = deltaTime
    self.fnameList = []

    # set arguments from the options
    self.verbose     = opts.verbose
    self.exttrig     = opts.followup_exttrig
    self.output_path = opts.output_path

    # default value for the injection time-window
    self.injectionWindow = 0.050

    # counter for the followups
    self.number = 0

    # for the estimated distances
    self.flow = opts.followup_flow
    self.spectrum = createSpectrum( self.flow, \
                                    sampleRate = 4096, \
                                    nPoints = 1048576)

    # retrieving the caches for the different stages
    self.triggerCache = {}
    for stage in self.stageLabels:
      pattern = stage
      self.triggerCache[stage] = cache.sieve(description=pattern)
      if self.opts.verbose:
        print "%d files found for stage %s" % (len(self.triggerCache[stage]), stage)

    pattern = "INJECTION_"
    self.injectionCache = cache.sieve(description=pattern).\
                          sieve(description=self.opts.followup_tag).\
                          sieve(ifos='HL')
  
    # generate a dictionary based on the event-ID in the exttrig case
    self.injections=dict()
    if self.exttrig:
      for file in self.injectionCache.pfnlist():
        injectionID = self.getInjectionID( file=file )    
        self.injections[injectionID] = SimInspiralUtils.\
                                       ReadSimInspiralFromFiles( [file], verbose=False )
    if self.verbose: print "parsing of cache files done..."

    # retrieve the time window used from one of the COIRE files
    if self.coireVetoMissedCache:
      coireFile = self.coireVetoMissedCache.pfnlist()[0]
      doc = utils.load_filename( coireFile, gz=(coireFile).endswith(".gz") )
      try:
        processParams = table.get_table(doc, lsctables.ProcessParamsTable.\
                                        tableName)
      except:
        raise "Error while reading process_params table from file", coireFile
    
      for tab in processParams:
        if tab.param=='--injection-window':
          self.injectionWindow = float(tab.value)/1000.0
      if self.verbose:
        print "Injection-window set to %.0f ms" % (1000*self.injectionWindow)
    else:
      self.injectionWindow = injectionWindow/1000.0

    # read the veto files
    self.readVetoFiles()


  # -----------------------------------------------------
  def getInjectionID( self, file=None, url=None ):
    """
    Extracting the injection ID from the filename, using
    the mechanism as used in lalapps_path2cache.
    (later: from a table in the file; not working now)
    @param file: filename from which the injection ID is extracted
    @param url:  url from which the injection ID is extracted
    """

    if file:
      path, file = os.path.split(file.strip())
      url = "file://localhost%s" % os.path.abspath(os.path.join(path, file))
      
    try:
      cache_entry = lal.CacheEntry.from_T050017(url)
    except ValueError, e:
      raise "Error while extracting injection ID from file ", filename

    injectionID = int( cache_entry.description.split('_')[-1] )
    return injectionID

  # -----------------------------------------------------
  def readVetoFiles( self ):
    """
    Reads the veto segments given by veto-files (if any)
    """
    self.vetodict = dict()

    # loop over the IFO names
    for ifoName in self.colors:

      self.vetodict[ifoName]=None
      # create the attribute name and check if it exists
      attributeName = 'followup_vetofile_'+ifoName.lower()
      if hasattr( self.opts, attributeName):

         # get the filename
         filename = getattr( self.opts, attributeName )

         if filename:
           #read the segment lists
           self.vetodict[ifoName] = segmentsUtils.fromsegwizard(open(filename))

  # -----------------------------------------------------
  def reset( self ):
    """
    Resets the counting number for the time-series plots generated.
    """
    self.number=0

  # -----------------------------------------------------
  def setTag( self, tag ):
    """
    Just sets the tag, called from plotinspmissed (why needed?)
    @param tag: tag for this followup
    """

    self.tag = tag

  # -----------------------------------------------------
  def print_inj( self, inj, injID ):
    """
    Print some useful informations to the screen.
    @param inj:   the current missed injection
    @param injID: the injection ID (used for exttrig only)
    """

    if self.exttrig:
      print "\nInjection details for injection %d with injID %s: " %\
            (self.number, injID)
    else:
      print "\nInjection details for injection %d:" % (self.number)    
      
    print "m1: %.2f  m2:%.2f  | end_time: %d.%d | "\
          "distance: %.2f  eff_dist_h: %.2f eff_dist_l: %.2f" % \
          ( inj.mass1, inj.mass2, self.end_time, \
	    self.end_time_ns, inj.distance, inj.eff_dist_h,\
	    inj.eff_dist_l )

  # ----------------------------------------------------
  def savePlot( self, stage ):
    """
    Saves the plots and store them in a seperate fnameList.
    @param stage: the stage this plot belongs to (e.g. INSPIRAL, THINCA,...)
    """
    fname = 'Images/'+self.opts.prefix + "_"+self.tag+"_map-"+\
            stage+"-"+str(self.number) +self.opts.suffix+'.png'
    fname_thumb = InspiralUtils.\
                  savefig_pylal( filename = self.output_path+fname,\
                                 doThumb = True, 
                                 dpi_thumb = self.opts.figure_resolution)
    
    self.fnameList.append( fname ) 
    return fname
 

  # -----------------------------------------------------  
  def findInjection( self, missedInj ):
    """
    Find the injection-ID corresponding to this particular
    missed injection.
    @param missedInj: the missed injection
    """
    
    # injID: the ID (a number)
    # groupInj: a list of SimInspiral tables
    for injID, groupInj in self.injections.iteritems():
      for inj in groupInj:
        if missedInj.geocent_end_time==inj.geocent_end_time and \
               missedInj.geocent_end_time_ns==inj.geocent_end_time_ns:
          return injID
          

    return None

  # -----------------------------------------------------  
  def getTimeTrigger( self, trig ):
    """
    This is a helper function to return a GPS time as one float number
    @param trig: a sngl_inspiral table entry
    """

    return float(trig.end_time) + float(trig.end_time_ns) * 1.0e-9
  

  # -----------------------------------------------------
  def getTimeSim(self, sim, ifo=None ):
    """
    This is a helper function to return a GPS time as one float number
    for a certain IFO. If no IFO is specified the injected geocentric
    time is returned.
    @param sim: a sim_inspiral table entry
    @param ifo: the IFO for which we want the sim time
    """
    
    time=0
    nano=0

    if not self.injection:
      time = self.end_time
      nano = self.end_time_ns
    elif not ifo:
      time = sim.geocent_end_time
      nano = sim.geocent_end_time_ns
    elif ifo:
      time = sim.get(ifo[0].lower()+'_end_time' )
      nano = sim.get(ifo[0].lower()+'_end_time_ns' )    

    return  float(time) + float(nano) * 1.0e-9

  # -----------------------------------------------------
  def highlightVeto( self, timeInjection, segLarge, ifoName, ylims  ):
    """
    Finds the intersection of the drawn triggers with the veto segments
    for this IFO
    """

    if not self.vetodict[ifoName]:
      return
    
    # find intersecting veto segments
    segVeto =  self.vetodict[ifoName] & segments.segmentlist( [segLarge] )

    # draw all intersecting segments
    for seg1 in segVeto:

      # after a tim-shift
      seg=seg1.shift( -timeInjection)
      vetox = [seg[0], seg[1], seg[1], seg[0], seg[0]]
      vetoy = [ylims[0], ylims[0], ylims[1], ylims[1], ylims[0]]
      fill ( vetox, vetoy, 'y', alpha=0.2)  


  # -----------------------------------------------------
  def isThereVeto( self, timeTrigger, ifoName ):
    
    flag = False
    vetoSegs = self.vetodict[ifoName]
    for seg in vetoSegs :
      if ( timeTrigger > seg[0] and timeTrigger < seg[1] ):
        flag = True

    return flag


  # -----------------------------------------------------
  def estimatedDistance( self, mass1, mass2, distTarget):

    snrActual = 8.0
    distActual = computeCandleDistance( 10.0, 10.0, \
                                        self.flow, self.spectrum, \
                                        snrActual)

    # rescale the SNR threshold
    snrTarget = distActual / distTarget*snrActual 
    distance = computeCandleDistance( mass1, mass2, \
                                      self.flow, self.spectrum,\
                                      snrTarget)
    return distance

  # -----------------------------------------------------
  def getExpectedSNR(self, triggerFiles, inj, number ):
    """
    Investigate template bank and returns exepcted horizon distance
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj: The current injection
    @param ifo: The current IFO
    @param number: The consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose: print "Processing TMPLTBANK triggers from files ",\
       triggerFiles      
    inspiralSumm, massInspiralSumm = InspiralUtils.\
                                     readHorizonDistanceFromSummValueTable(\
                                     triggerFiles, self.verbose)

    # get different informations on the missed injection
    injMass = [inj.mass1, inj.mass2]
    timeInjection = self.getTimeSim( inj )    
    totalMass  = inj.mass1 + inj.mass2
    eta = inj.mass1 * inj.mass2 / totalMass / totalMass

    # Given the masses (and therefore eta), the expected horizon distance(SNR) has to be scaled by this factor
    factor = sqrt(4 * eta)

    output = {}
    # for each ifo that has been found
    for ifo in massInspiralSumm.keys():     
      # loop over the summ value table as long as horizon is not set
      output[ifo] = [] 
      horizon = 0 
      for massNum in range(len(massInspiralSumm[ifo])):
        # looking at the masses and horizon distance (and threshold)
        if horizon > 0:
          break
        for this, horizon  in zip(massInspiralSumm[ifo][massNum].\
                                  getColumnByName('comment'),\
            massInspiralSumm[ifo][massNum].getColumnByName('value').asarray()):
          masses = this.split('_')
          threshold = float(masses[2])
          readTotalMass = float(masses[0])+float(masses[1])
          # check if the current total Mass is greater tan the requested total Mass
          # if so, then keep this horizon and break 
          if readTotalMass > totalMass:
            startTimeSec = \
              massInspiralSumm[ifo][massNum].getColumnByName('start_time').asarray()
            #sanity check 
            if len(startTimeSec)>1: 
              print >> sys.stderr, 'Warning in fct. expectedHorizonDistance: '\
                    'More than 1 file found at particular GPS time. Using the first file.' 
            if startTimeSec[0] > self.end_time:
              text= """the start time of the template bank must be less than 
		the end_time of the injection. We found startTime of the bank to be %s and
		 the geocent end time of the injection to be %s""",startTimeSec[0], self.end_time
              raise ValueError,  text
	    if self.injection:
              output[ifo] = horizon * factor * threshold / float(getattr(inj, 'eff_dist_'+ifo[0].lower() ))
	    else:
	      output[ifo] = horizon * factor * threshold / float(getattr(inj, 'eff_distance' ))
            break
          #'otherwise, reset horizon distance
          else: horizon = 0           

    return output

  # -----------------------------------------------------
  def putText( self, text ):
    """
    Puts some text into an otherwise empty plot.
    @param text: text to put in the empty plot
    """
    newText = ''
    for i in range( int(len(text)/60.0)+1):
      newText+=text[60*i:60*i+60]+'\n'
    figtext(0.15,0.15, newText)
 
  # -----------------------------------------------------
  def investigateTimeseries(self, triggerFiles, inj,  ifoName, stage, number,\
                            yval = 'snr'):
    """
    Investigate inspiral triggers and create a time-series
    of the SNRs around the injected time
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj:          the current missed injection
    @param ifoName:      the IFO for which the plot is made 
    @param stage:        the name of the stage (FIRST, SECOND)
    @param number:        the consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose: print "Processing INSPIRAL triggers from files ", triggerFiles   
    snglTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
      triggerFiles , verbose=False)

    # create a figure
    fig=figure()
    foundSet = set()
    trigger_details = {}
    loudest_details = {}
    
    if snglTriggers is None:
      # put message on the plot instead
      #print "NO SNGL TRIGGERS HERE"
      self.putText( 'No sngl_inspiral triggers in %s' % str(triggerFiles))

    else:
      # selection segment
      timeInjection = self.getTimeSim( inj )
      segSmall =  segments.segment( timeInjection-self.injectionWindow, \
                                    timeInjection+self.injectionWindow )
      segLarge =  segments.segment( timeInjection-self.deltaTime, \
                                    timeInjection+self.deltaTime )

      # create coincidences for THINCA stage
      coincTriggers = None
      if 'THINCA' in stage:
        coincTriggers = CoincInspiralUtils.coincInspiralTable( snglTriggers, \
                      CoincInspiralUtils.coincStatistic("snr") )
        selectedCoincs = coincTriggers.vetoed( segSmall )
      
      # loop over the IFOs (although this is a plot for IFO 'ifoName')
      for ifo in self.colors.keys():

        # get the singles for this ifo
        snglInspiral = snglTriggers.ifocut(ifo)

        # select a range of triggers
        selectedLarge = snglInspiral.vetoed( segLarge )
        timeLarge = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedLarge ]

        selectedSmall = snglInspiral.vetoed( segSmall )
        timeSmall = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedSmall ]

        # use the set of selected coincident triggers in the THINCA stages
        if coincTriggers:
#          selectedSmall = selectedCoincs.getsngls(ifo)
          selectedSmall = selectedCoincs.cluster(2* self.injectionWindow).getsngls(ifo)
          timeSmall = [ self.getTimeTrigger( sel )-timeInjection \
                        for sel in selectedSmall ]
          
        # skip if no triggers in the large time window
        if len(timeLarge)==0:
          self.putText( 'No triggers/coincidences found within time window')
          continue

        # add IFO to this set; the injection is found for this IFO and stage
        if len(timeSmall)>0:
          foundSet.add(ifo)
          
          # record some details from the 50ms triggers
          trigger_details[ifo] = {}
          trigger_details[ifo]["snr"] = selectedSmall.get_column('snr')
          trigger_details[ifo]["chisq"] = selectedSmall.get_column('chisq')
          trigger_details[ifo]["eff_dist"] = selectedSmall.get_column('eff_distance')
          trigger_details[ifo]["mchirp"] = selectedSmall.get_column('mchirp')

          # record details of the loudest trigger
          loudest_details[ifo] = {}
          loudest = selectedSmall[selectedSmall.get_column('snr').argmax()]
          loudest_details[ifo]["snr"] = loudest.snr
          loudest_details[ifo]["mchirp"] = loudest.mchirp
          loudest_details[ifo]["eff_dist"] = loudest.eff_distance
          loudest_details[ifo]["chisq"] = loudest.chisq
          loudest_details[ifo]["timeTrigger"] = self.getTimeTrigger( loudest )

          timeTrigger = self.getTimeTrigger( loudest )
          vetoSegs = self.vetodict[ifoName]
          
        # plot the triggers
        plot( timeLarge, selectedLarge.get_column(yval),\
              self.colors[ifo]+'o', label="_nolegend_")
        plot( timeSmall, selectedSmall.get_column(yval), \
              self.colors[ifo]+'s', label=ifo)

      # draw the injection times and other stuff
      ylims=axes().get_ylim()
      plot( [0,0], ylims, 'g--', label="_nolegend_")
      plot( [-self.injectionWindow, -self.injectionWindow], ylims, 'c:',\
            label="_nolegend_")
      plot( [+self.injectionWindow, +self.injectionWindow], ylims, 'c:',\
            label="_nolegend_")

      self.highlightVeto( timeInjection, segLarge, ifoName, ylims  )

      # save the plot
      grid(True)
      legend()

    ylims=axes().get_ylim()
    axis([-self.deltaTime, +self.deltaTime, ylims[0], ylims[1]])
    xlabel('time [s]')
    ylabel(yval)
    title(stage+'_'+str(self.number))    
    fname = self.savePlot( stage )
    close(fig)

    result = {'filename':fname, 'foundset':foundSet, 'trigger_details':trigger_details, 'loudest_details':loudest_details}
    return result

  # -----------------------------------------------------
  def investigateTrigbank(self, triggerFiles, inj, ifoName,  stage, number ):
    """
    Investigate the trigbank
    @param triggerFiles: List of files containing the trigbank
    @param inj: The current injection
    @param ifo: The current IFO    
    @param stage: The name of the stage (FIRST, SECOND)
    @param number: The consecutive number for this inspiral followup
    """
    
    if self.injection:
      injMass = [inj.mass1, inj.mass2]
    elif self.coincInspiral:
      injMass = [self.Sngls[ifoName].mass1,self.Sngls[ifoName].mass2]

    # read the trigbank
    if self.verbose: print "Processing TRIGBANK triggers from files ", triggerFiles
    snglTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
      triggerFiles, verbose=False )

    fig=figure()
    selected = dict()
    
    if snglTriggers is None:
      self.putText( 'no sngl inspiral triggers found in %s' % str(triggerFiles))

    else:
      for ifo in self.colors.keys():

        # get the singles and plot the tmpltbank
        snglInspiral = snglTriggers.ifocut(ifo)
        plot( snglInspiral.get_column('mass1'), \
              snglInspiral.get_column('mass2'), self.colors[ifo]+'x', \
              label=ifo)

    # plot the missed trigger and save the plot
    xm = injMass[0]
    ym = injMass[1]
    plot( [xm], [ym],  'ko', ms=10.0, mfc='None', mec='k', mew=3, \
          label="_nolegend_")
    grid(True)
    legend()
    xlabel('mass1')
    ylabel('mass2')
    title(stage+'_'+str(self.number))    
    fname = self.savePlot( stage )
    close(fig)

    result = {'filename':fname, 'foundset':None}
    return result

  # -----------------------------------------------------    
  def selectCategory(self, triggerFiles, category):
    """
    Return a trigger list that contains only files for the choosen category.
    @param triggerList : a list of file names
    @param category: a category number
    @return: a sub list of filename corresponding to the category requested
    """
    
    if category==1:
      # for now (June 2008, TC) CAT_1 is not part of the filename,
      # so choose only files without 'CAT' in the filename
      newList = [file for file in triggerFiles if 'CAT' not in file]
                     
    else:
      cat = 'CAT_'+str(category)
      newList = [file for file in triggerFiles if cat in file]      
          
    return newList


  # -----------------------------------------------------
  def followup(self, inj, selectIFO, description=None,type='SimInspiral' ):
    """
    Do the followup procedure for the missed injection 'inj'
    and create the several time-series for INSPIRAL, THINCA and
    the trigbank.
    The return value is the name of the created html file
    @param inj: sim_inspiral table of the missed injection
    @param selectIFO: the IFO that is investigated
    @param description : a description to select files from the cache file (followup_tag).
    """
    
    def fillTable(page, contents ):
      """
      Making life easier...
      """
      page.add('<tr>')
      for content in contents:
        page.add('<td>')
        page.add( str(content) )
        page.add('</td>')
      page.add('</tr>')


    def fillTableMore(page, contents, moreContents ):
      """
      Making life easier...
      """
      page.add('<tr>')
      for content in contents:
        page.add('<td>')
        page.add( str(content) )
        page.add('</td>')
      page.add('</tr>')

   
    # get the injections that are missed
    injID = 0
    if self.exttrig and type == 'SimInspiral':
      injID = self.findInjection( inj )

    if type == 'SimInspiral':
      self.injection = True
      self.end_time = inj.geocent_end_time
      self.end_time_ns = inj.geocent_end_time_ns
    elif type == 'coincInspiral':
      self.coincInspiral = True
      self.h1 = False
      self.h2 = False
      self.l1 = False
      self.v1 = False
      self.Sngls = {}
      self.ifos = []
      self.injection = False
      if 'H1' in inj.get_ifos()[1]:
        self.ifos.append('H1')
        self.Sngls['H1'] = inj.H1
      if 'H2' in inj.get_ifos()[1]:
        self.ifos.append('H2')
        self.Sngls['H2'] = inj.H2
      if 'L1' in inj.get_ifos()[1]:
        self.ifos.append('L1')
        self.Sngls['L1'] = inj.L1
      if 'V1' in inj.get_ifos()[1]:
        self.ifos.append('V1')
        self.Sngls['V1'] = inj.V1
      self.end_time = self.Sngls[selectIFO].end_time
      self.end_time_ns = self.Sngls[selectIFO].end_time_ns
    else:
      print >> sys.stderr, 'Type not recognized, please enter "SimInspiral'
      print >> sys.stderr, 'or "coincInspiral" in the type field'
    # increase internal number:
    self.number+=1

    ## create the web-page
    page = markup.page()
    page.h1("Followup missed injection #"+str(self.number)+" in "+selectIFO )
    page.hr()

    # Check if this is a spin Taylor injection and calculate spin values

    if self.injection:
      if 'SpinTaylor' in inj.waveform:
        chi1 = sqrt(inj.spin1x*inj.spin1x \
                    +inj.spin1y*inj.spin1y+inj.spin1z*inj.spin1z)
        chi2 = sqrt(inj.spin2x*inj.spin2x \
                    +inj.spin2y*inj.spin2y+inj.spin2z*inj.spin2z)
    
    # add a table
    page.add('<table border="3" ><tr><td>')
    page.add('<table border="2" >')          
    fillTable( page, ['<b>parameter','<b>value'] )
    fillTable( page, ['Number', self.number] )
    if self.exttrig:
      fillTable( page, ['inj ID', injID] )
    if self.injection:
      fillTable( page, ['mass1', '%.2f'% inj.mass1] )
      fillTable( page, ['mass2', '%.2f'%inj.mass2] )
      fillTable( page, ['mtotal', '%.2f' % (inj.mass1+inj.mass2)] )
      fillTable( page, ['mchirp', '%.2f' % (inj.mchirp)] )
      if 'SpinTaylor' in inj.waveform:
        fillTable( page, ['chi1', '%.4f' % chi1] )
        fillTable( page, ['chi2', '%.4f' % chi2] )
        fillTable( page, ['spin1x', '%.4f' % inj.spin1x] )
        fillTable( page, ['spin1y', '%.4f' % inj.spin1y] )
        fillTable( page, ['spin1z', '%.4f' % inj.spin1z] )
        fillTable( page, ['spin2x', '%.4f' % inj.spin2x] )
        fillTable( page, ['spin2y', '%.4f' % inj.spin2y] )
        fillTable( page, ['spin2z', '%.4f' % inj.spin2z] )
        fillTable( page, ['theta0', '%.4f' % inj.theta0] )
        fillTable( page, ['phi0', '%.4f' % inj.phi0] )
      fillTable( page, ['end_time', self.end_time] )
      fillTable( page, ['end_time_ns', self.end_time_ns] )    
      fillTable( page, ['distance', '%.1f' % inj.distance] )
      fillTable( page, ['eff_dist_h','%.1f' %  inj.eff_dist_h] )
      fillTable( page, ['eff_dist_l','%.1f' %  inj.eff_dist_l] )
      fillTable( page, ['eff_dist_v','%.1f' %  inj.eff_dist_v] )
      fillTable( page, ['eff_dist_g','%.1f' %  inj.eff_dist_g] )  
      fillTable( page, ['playground','%s' %  InspiralUtils.isPlayground(inj)] )
    else:
      if type == 'coincInspiral':
        for ifo in self.ifos:
          fillTable( page, ['mass1 ' + ifo, '%.2f'% self.Sngls[ifo].mass1] )
          fillTable( page, ['mass2 ' + ifo, '%.2f'% self.Sngls[ifo].mass2] )
          fillTable( page, ['mchirp ' + ifo, '%.2f'% self.Sngls[ifo].mchirp] )
        for ifo in self.ifos:
          fillTable( page, ['end time ' + ifo, \
                            self.Sngls[ifo].end_time] )
          fillTable( page, ['end time (ns) ' + ifo, \
                            self.Sngls[ifo].end_time_ns] )
        for ifo in self.ifos:
          fillTable( page, ['Effective distance ' + ifo, \
                            '%.2f'% self.Sngls[ifo].eff_distance] )
        
    page.add('</table></td>')
    
    # print infos to screen as well
    if self.injection:
      if self.opts.verbose: self.print_inj( inj,  injID)

    # retrieve other files for this missed injection
    investDict = {}
    for stage, cache in self.triggerCache.iteritems():

      trigCache = lal.Cache()
      for c in cache:
        if self.end_time in c.segment:
          if self.exttrig and 'TMPLTBANK' not in stage:
            if self.getInjectionID( url = c.url )==injID:
              trigCache.append( c )
          else:
            trigCache.append( c )            

      # create a filelist
      fileList = trigCache.sieve(description=description).pfnlist()
      if 'TRIGBANK' in stage:
        fileList=trigCache.sieve(description=description).sieve(ifos=selectIFO).pfnlist()        
        
      # check if the pfnlist is empty. `
      if len(fileList)==0:
        print >>sys.stderr, stage,injId,self.end_time
        print >>sys.stderr, "Error: No files found for stage %s in the "\
              "cache for ID %s and time %d; probably mismatch of a "\
              "pattern in the options. " % \
              ( stage, injID, self.end_time)        
        continue

      # now create several plots
      if 'TMPLTBANK' in stage:
        investDict[stage]= self.investigateTrigbank( fileList, inj, selectIFO, stage, self.number )
	if self.injection:
          investDict[stage]['horizon'] = self.getExpectedSNR( fileList, inj, self.number )
      elif 'TRIGBANK' in stage:
        investDict[stage] = self.investigateTrigbank( fileList, inj, selectIFO, stage, self.number )
      elif 'THINCA_SECOND' in stage:

        # loop over the four categories
        for cat in [1,2,3,4]:
          selectList=self.selectCategory( fileList, cat)
          if len(selectList)==0: continue
          modstage=stage+'_CAT_'+str(cat)
          investDict[modstage] = self.investigateTimeseries( selectList, inj, selectIFO, modstage, self.number )
          investDict[modstage+'_mchirp'] = self.investigateTimeseries(\
              selectList, inj, selectIFO, modstage + '_mchirp', \
              self.number,yval='mchirp' )
      else:
        investDict[stage]=self.investigateTimeseries( fileList, inj, selectIFO, stage, self.number)
        investDict[stage+ '_mchirp']=self.investigateTimeseries( fileList,\
	    inj, selectIFO, stage + '_mchirp', self.number,yval='mchirp')
          
    ## print out the result for this particular injection
    page.add('<td><table border="2" >')
    fillTable( page, ['<b>step','<b>F/M', '<b>Rec. SNR', '<b>Rec. mchirp', '<b>Rec. eff_dist', '<b>Rec. chisq', '<b>Veto ON/OFF'] )
    
    for stage in self.orderLabels:
      if investDict.has_key(stage):
        result = investDict[stage]

        # create statements on the followup page
        if "TRIGBANK" in stage:
          pass 
      
        elif "TMPLTBANK" in stage:
	  if self.injection:
            for ifo in result['horizon'].keys():
              text = "Expected SNR for %s is " % ifo
              this_snr = "%.2f" % result['horizon'][ifo]
              fillTable( page, [ text,  this_snr]) 

        elif result['foundset']:
          foundIFOs=''
          if "INSPIRAL" in stage or "THINCA" in stage:
            foundIFOs=' in '
            recSnr=''
            recChirpMass=''
            recEffDist=''
            recChisq=''
            loudestSnr=''
            loudestChirpMass=''
            loudestEffDist=''
            loudestChisq=''
            vetoOnOff=''
            for ifo in result['foundset']:
              foundIFOs+=ifo+' '
              # Pars of the 50ms recovered triggers
              recSnr+=ifo+str(result['trigger_details'][ifo]['snr'])+'<br>'
              recChirpMass+=ifo+str(result['trigger_details'][ifo]['mchirp'])+'<br>'
              recEffDist+=ifo+str(result['trigger_details'][ifo]['eff_dist'])+'<br>'
              recChisq+=ifo+str(result['trigger_details'][ifo]['chisq'])+'<br>'
              # Pars of the loudest trigger
              loudestSnr+=ifo+': '+str(result['loudest_details'][ifo]['snr'])+'<br>'
              loudestChirpMass+=ifo+': '+str(result['loudest_details'][ifo]['mchirp'])+'<br>'
              loudestEffDist+=ifo+': '+str(result['loudest_details'][ifo]['eff_dist'])+'<br>'
              loudestChisq+=ifo+': '+str(result['loudest_details'][ifo]['chisq'])+'<br>'
              # Whether some of the ifo times is vetoed
              timeTrigger = float(result['loudest_details'][ifo]['timeTrigger'])
	      if (self.vetodict[ifo]):
                Veto = self.isThereVeto ( timeTrigger, ifo )
                if (Veto):
                  onOff = 'ON'
                  vetoOnOff+=ifo+': '+onOff+'<br>'
                else:
                  onOff = 'OFF'
                  vetoOnOff+=ifo+': '+onOff+'<br>'
	      else: 
		  vetoOnOff+=ifo+': No info<br>'
          fillTable( page, [ stage,  'FOUND'+foundIFOs, 'loudest<br>'+loudestSnr, 'loudest<br>'+loudestChirpMass, 'loudest<br>'+loudestEffDist, 'loudest<br>'+loudestChisq, vetoOnOff])

        else:
          fillTable( page, [ stage,  '<font color="red">MISSED'])      
    page.add('</table>')
    page.add('</td></tr></table><br><br>')


    ## add the pictures to the webpage
    for stage in self.orderLabels:
      if investDict.has_key(stage):
        result = investDict[stage]
      
        if stage!="TMPLTBANK":
          fname = result['filename']
          page.a(extra.img(src=[fname], width=400, \
                           alt=fname, border="2"), title=fname, href=[ fname ])
      if investDict.has_key(stage + '_mchirp'):
        result = investDict[stage + '_mchirp']
        if stage!="TMPLTBANK":
          fname = result['filename']
          page.a(extra.img(src=[fname], width=400, \
                           alt=fname, border="2"), title=fname, href=[ fname ])
        
    # and write the html file
    htmlfilename = self.opts.prefix + "_"+selectIFO+"_followup_"+str(self.number) +\
                         self.opts.suffix+'.html'
    file = open(self.opts.output_path+htmlfilename,'w')      
    file.write(page(False))
    file.close()

    # store html file in fnameList
    self.fnameList.append(htmlfilename)
    
    # supply the output
    return htmlfilename


# -----------------------------------------------------
def createSpectrum( flow, sampleRate = 4096, nPoints = 1048576):

  def calcPSD( f ):
    f0=150.0
    a = pow(4.49*f/f0, -56.0)
    b = 0.16*pow(f/f0, -4.52)
    c = 0.52
    d = 0.32*pow(f/f0,2)
    return (a+b+c+d)*9.0e-46
  
  chanDeltaT = 1.0/float(sampleRate)
  deltaF =  1.0/( float(nPoints)*chanDeltaT )
  
  spectrum = []
  for k in range( nPoints/2+1):
    f = k*deltaF
    if f<flow:
      spectrum.append(0)
    else:
      spectrum.append( calcPSD(f) )
    
  return spectrum
      
# -----------------------------------------------------
def computeCandleDistance( candleM1, candleM2, fLow, spectrum = None, \
                           snr=8.0, sampleRate = 4096, nPoints = 1048576):
  """
  Computes the candle distance as computed in inspiralutils.c.
  This code has been tested on 21 Feb 2008 and the derivation between
  these calculated values with the LAL code is less than 2 Percent.
  """
  # TestCode: testEstimatedDistance.png
  

  LAL_MTSUN_SI = 4.92549095e-6
  chanDeltaT = 1.0/float(sampleRate)
  negativeSevenOverThree = -7.0/3.0
  totalMass = candleM1 + candleM2
  mu = candleM1 * candleM2 / totalMass
  distNorm = 9.5708317e-20 # = 2.0 * LAL_MRSUN_SI / (1.0e6 * LAL_PC_SI )
  a = sqrt( (5.0 * mu) / 96.0 ) * \
      pow( totalMass / ( pi*pi ), 1.0/3.0 ) *\
      pow( LAL_MTSUN_SI / chanDeltaT, -1.0/6.0 )
  sigmaSq = 4.0 * ( chanDeltaT / float(nPoints) ) * \
            distNorm * distNorm * a * a
  
  fmax = 1.0 / (6.0 * sqrt(6.0) * pi * totalMass * LAL_MTSUN_SI)
  deltaF =  1.0/( float(nPoints)*chanDeltaT )

  cut = int( fLow/deltaF )
  kmax = min( int( fmax/deltaF ), len(spectrum) )

  sigmaSqSum = 0
  for k in range( cut, kmax):
    sigmaSqSum += pow( float(k)/float(nPoints),  negativeSevenOverThree ) /\
                  spectrum[k]

  sigmaSq *= sigmaSqSum
  distance = sqrt( sigmaSq ) / snr

  return distance

