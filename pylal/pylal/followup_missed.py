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

import os, sys, exceptions
from optparse import *

import matplotlib
matplotlib.use('Agg')
from pylab import *
try:
  set
except NameError:
  from sets import Set as set

from pylal import SnglInspiralUtils
from pylal import SimInspiralUtils
from pylal import CoincInspiralUtils
from pylal import InspiralUtils
from glue import lal
from glue import segments

from pylal import webUtils
from pylal import readMeta
from glue import markup
from glue.markup import oneliner as extra
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import segments
from glue import segmentsUtils

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
  def __init__(self, cache, coireVetoMissedCache, opts, deltaTime=20):
    """
    Initialize this class and sets up all the cache files.
    @param cache: the cache of all files
    @param coireVetoMissedCache: the cache that contains the COIRE-files
                                  of the missed injections
    @param opts: The 'opts' from the main code
    @param deltaTime: The time window that is used in the time-series plots (in seconds)
    """
    from glue.ligolw import table # must do this within here for unknown reason; error instead
    rcParams.update({'text.usetex': False})

    self.colors = {'H1':'r','H2':'b','L1':'g','V1':'m','G1':'c'}
    self.stageLabels = ['TMPLTBANK', 'INSPIRAL_FIRST','THINCA_FIRST','TRIGBANK', 'INSPIRAL_SECOND', \
	'THINCA_SECOND', 'THINCA_SECOND_CAT_2', 'THINCA_SECOND_CAT_3', 'THINCA_SECOND_CAT_4']

    # setting all the parameters
    self.cache                = cache
    self.coireVetoMissedCache = coireVetoMissedCache
    self.opts                 = opts
    self.deltaTime = deltaTime
    self.fnameDict = dict()

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

    # getting all the caches
    self.stageLabels = ['TMPLTBANK', 'INSPIRAL_FIRST','THINCA_FIRST',\
                        'TRIGBANK', 'INSPIRAL_SECOND', 'THINCA_SECOND']
    self.triggerCache = {}
    for stage in self.stageLabels:
      pattern = stage
      if self.opts.user_tag:
        pattern+='*_'+self.opts.user_tag
      self.triggerCache[stage] = cache.sieve(description=pattern)

    pattern = "INJECTION_"
    if self.opts.user_tag:
      pattern+='*'+self.opts.user_tag
    self.injectionCache = cache.sieve(description=pattern).sieve(ifos='HL')
  
    # generate a dictionary based on the event-ID
    # in the exttrig case
    self.injections=dict()
    if self.exttrig:
      for file in self.injectionCache.pfnlist():
        basename = os.path.basename(file)
        
        injectionString =  basename.split(self.opts.user_tag)[1].\
                          split('_')[1].split('-')[0]
        injectionID =  int(float(injectionString))
        self.injections[injectionID] = SimInspiralUtils.\
                                       ReadSimInspiralFromFiles( [file] )
        p = readMeta.metaDataTable( [ file ] , 'process_params')
    if self.verbose:
      print "parsing of cache files done..."

    # retrieve the time window used from one of the COIRE files
    coireFile = self.coireVetoMissedCache.pfnlist()[0]
    doc = utils.load_filename( coireFile, gz=(coireFile).endswith(".gz") )
    try:
      processParams = table.get_table(doc, lsctables.ProcessParamsTable.\
                                      tableName)
    except:
      print >>sys.stderr,"Error while reading process_params table from file ",\
            coireFile
      print >>sys.stderr, "Probably files does not exist"
      sys.exit(1)

    for table in processParams:
      if table.param=='--injection-window':
        self.injectionWindow = float(table.value)/1000.0
    if self.verbose:
      print "Injection-window set to %.0f ms" % (1000*self.injectionWindow)


    self.readVetoFiles()
    

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
    Resets the counting number
    """
    self.number=0

  # -----------------------------------------------------
  def setTag( self, tag ):
    """
    Just sets the tag...
    """

    self.tag = tag

  # -----------------------------------------------------
  def print_inj( self, inj, injID ):
    """
    Print some useful informations to the screen if the verbose
    flag is set.
    @param inj: the current missed injection
    @param injID: the injection ID (useful for exttrig only)
    """

    if self.exttrig:
      print "\nInjection details for injection %d with injID %s: " % (self.number, injID)
    else:
      if self.verbose is True :
        print "\nInjection details for injection %d:" % (self.number)    
    if self.verbose is True:
      print "m1: %.2f  m2:%.2f  | end_time: %d.%d | "\
        "distance: %.2f  eff_dist_h: %.2f eff_dist_l: %.2f" % \
        ( inj.mass1, inj.mass2, inj.geocent_end_time, inj.geocent_end_time_ns,\
          inj.distance, inj.eff_dist_h, inj.eff_dist_l )

  # ----------------------------------------------------
  def savePlot( self, stage ):
    """
    Saves the plots and store them in a seperate fnameList
    """
    fname = 'Images/'+self.opts.prefix + "_"+self.tag+"_map-"+\
            stage+"-"+str(self.number) +self.opts.suffix+'.png'
    fname_thumb = InspiralUtils.\
                  savefig_pylal( filename = self.output_path+fname,\
                                 doThumb = True, 
                                 dpi_thumb = self.opts.figure_resolution)

    self.fnameDict[stage]=fname

  # -----------------------------------------------------  
  def findInjection( self, missedInj ):
    """
    Find the injection-ID corresponding to this particular
    missed injection 'missedInj'
    @param missedInj: This missed injection
    """
    
    # injID: one ID (a number)
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
    Just to make life easier...
    @param trig: A sngl_inspiral table
    """

    return float(trig.end_time) + float(trig.end_time_ns) * 1.0e-9
  

  # -----------------------------------------------------
  def getTimeSim(self, sim, ifo=None ):
    """
    Just to make life easier...
    @param sim: A sim_inspiral table
    @param ifo: The IFO
    """
    
    time=0
    nano=0

    if not ifo:
      time = sim.geocent_end_time
      nano = sim.geocent_end_time_ns
    if ifo:
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
  def expectedHorizonDistance(self, triggerFiles, inj,  ifoName, number ):
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
     
    injMass = [inj.mass1, inj.mass2]
    # selection segment
    timeInjection = self.getTimeSim( inj )
    segSmall =  segments.segment( timeInjection-self.injectionWindow, \
                                  timeInjection+self.injectionWindow )
    segLarge =  segments.segment( timeInjection-self.deltaTime, \
                                  timeInjection+self.deltaTime )
    foundAny = set()
    totalMass  = inj.mass1 + inj.mass2
    eta = inj.mass1 * inj.mass2 / totalMass / totalMass

    # Given the masses (and therefore eta), the expected horizon distance(SNR) has to be scaled by this factor
    factor = math.sqrt(4 * eta)

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
              print >> sys.stderr, 'Warning in fct. expectedHorizonDistance: More than 1 file found at particular GPS time. Using the first file.' 
            if startTimeSec[0] > inj.geocent_end_time:
              text= """the start time of the template bank must be less than 
		the end_time of the injection. We found startTime of the bank to be %s and
		 the geocent end time of the injection to be %s""",startTimeSec[0], inj.geocent_end_time
              raise ValueError,  text
            output[ifo] = horizon * factor * threshold / float(getattr(inj, 'eff_dist_'+ifo[0].lower() ))
            break
          #'otherwise, reset horizon distance
          else: horizon = 0 
            


    return output

  # -----------------------------------------------------
  def investigateInspiral(self, triggerFiles, inj,  ifoName, stage, number ):
    """
    Investigate inspiral triggers and create a time-series
    of the SNRs around the injected time
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj: The current injection
    @param ifo: The current IFO
    @param stage: The name of the stage (FIRST, SECOND)
    @param number: The consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose: print "Processing INSPIRAL triggers from files ", triggerFiles   
    snglTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
      triggerFiles , verbose=self.verbose)

    if snglTriggers is None:
      figtext(0,0,'no sngl_inspiral table in %s' % str(triggerFiles))

    else:
      # selection segment
      timeInjection = self.getTimeSim( inj )
      segSmall =  segments.segment( timeInjection-self.injectionWindow, \
                                    timeInjection+self.injectionWindow )
      segLarge =  segments.segment( timeInjection-self.deltaTime, \
                                    timeInjection+self.deltaTime )

      foundAny = set()
      fig=figure()
      for ifo in self.colors.keys():
        # get the singles 
        snglInspiral = snglTriggers.ifocut(ifo)
      
        # select a range of triggers
        selectedLarge = snglInspiral.vetoed( segLarge )
        timeLarge = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedLarge ]

        selectedSmall = snglInspiral.vetoed( segSmall )
        timeSmall = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedSmall ]

        if len(timeLarge)==0:
          continue

        if len(timeSmall)>0:
          foundAny.add(ifo)


        # plot the triggers
        plot( timeLarge, selectedLarge.get_column('snr'), self.colors[ifo]+'o',\
              label="_nolegend_")
        plot( timeSmall, selectedSmall.get_column('snr'), self.colors[ifo]+'s', \
            label=ifo)

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
      
    axis([-self.deltaTime, +self.deltaTime, ylims[0], ylims[1]])
    xlabel('time [s]')
    ylabel('SNR')
    title(stage+'_'+str(self.number))    
    self.savePlot( stage )
    close(fig)

    return foundAny

  # -----------------------------------------------------
  def investigateThinca(self, triggerFiles, inj, ifoName, stage, number ):
    """
    Investigate coincdent triggers and create a time-series
    of the SNRs around the injected time
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj: The current injection
    @param ifo: The current IFO    
    @param stage: The name of the stage (FIRST, SECOND)
    @param number: The consecutive number for this inspiral followup
    """
    
    # read and construct the coincident events
    if self.verbose: print "Processing THINCA triggers from files "
    snglInspiral = SnglInspiralUtils.ReadSnglInspiralFromFiles( triggerFiles, \
	mangle_event_id=True, verbose=self.verbose )
    coincInspiral = CoincInspiralUtils.coincInspiralTable( snglInspiral, \
	CoincInspiralUtils.coincStatistic("snr") )

    if snglInspiral is None:
      figtext(0,0,'no sngl_inspiral table in %s' % str(triggerFiles))
      
    else:

      # selection segment
      timeInjection = self.getTimeSim( inj )
      segSmall =  segments.segment( timeInjection-self.injectionWindow, \
                                    timeInjection+self.injectionWindow )
      segLarge =  segments.segment( timeInjection-self.deltaTime, \
                                    timeInjection+self.deltaTime )

      plotAny = False
      foundAny= False
      fig=figure()
      selected = dict()  
      for ifo in self.colors.keys():

        # get the singles 
        snglInspiral = coincInspiral.getsngls(ifo)    
        
        # select a range of triggers
        selectedLarge = snglInspiral.vetoed( segLarge )
        timeLarge = [ self.getTimeTrigger( sel )-timeInjection\
                      for sel in selectedLarge ]

        selectedSmall = snglInspiral.vetoed( segSmall )
        timeSmall = [ self.getTimeTrigger( sel )-timeInjection\
                      for sel in selectedSmall ]

        if len(timeLarge)==0:
          continue

        if len(timeSmall)>0:
          foundAny=True

        # plot the triggers
        plotAny = True
        plot( timeLarge, selectedLarge.get_column('snr'), self.colors[ifo]+'o',\
              label="_nolegend_")
        plot( timeSmall, selectedSmall.get_column('snr'), self.colors[ifo]+'s', \
            label=ifo)


      # draw the injection times and other stuff
      if plotAny:
        ylims=axes().get_ylim()
      else:
        ylims=[0,1]
      plot( [0,0], ylims, 'g--', label="_nolegend_")
      plot( [-self.injectionWindow, -self.injectionWindow], ylims, 'c:',\
            label="_nolegend_")
      plot( [+self.injectionWindow, +self.injectionWindow], ylims, 'c:',\
            label="_nolegend_")

      # save the plot
      grid(True)
      if plotAny:
        legend()
      axis([-self.deltaTime, +self.deltaTime, ylims[0], ylims[1]])
      
    xlabel('time [s]')
    ylabel('SNR')
    title(stage+'_'+str(self.number))
    self.savePlot( stage )
    close(fig)

    return foundAny
  
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

    injMass = [inj.mass1, inj.mass2]

    # read the trigbank
    if self.verbose: print "Processing TRIGBANK triggers from files "
    snglTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
      triggerFiles, verbose=self.verbose )


    fig=figure()
    selected = dict()  
    for ifo in self.colors.keys():

      # get the singles 
      snglInspiral = snglTriggers.ifocut(ifo)

      # plot the templatebank
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
    self.savePlot( stage )
    close(fig)


  # -----------------------------------------------------
  def followup(self, inj, ifo ,description=None):
    """
    Do the followup procedure for the missed injection 'inj'
    and create the several time-series for INSPIRAL, THINCA and
    the trigbank.
    The return value is the name of the created html file
    @param inj: sim_inspiral table of the missed injection
    @param ifo : The current IFO
    @param ifo : a description to select files in a cache file. It could be the tag of the injection run.
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

   
    # get the injections that were missed
    injID = 0
    if self.exttrig:
      injID = self.findInjection( inj )

    # increase number:
    self.number+=1

    ## create the web-page
    page = markup.page()
    page.h1("Followup missed injection #"+str(self.number)+" in "+ifo )
    page.hr()
    
    # add a table
    page.add('<table border="3" ><tr><td>')
    page.add('<table border="2" >')          
    fillTable( page, ['<b>parameter','<b>value'] )
    fillTable( page, ['Number', self.number] )
    fillTable( page, ['inj ID', injID] )
    fillTable( page, ['mass1', '%.2f'% inj.mass1] )
    fillTable( page, ['mass2', '%.2f'%inj.mass2] )
    fillTable( page, ['mtotal', '%.2f' % (inj.mass1+inj.mass2)] )
    fillTable( page, ['mchirp', '%.2f' % (inj.mchirp)] )
    fillTable( page, ['end_time', inj.geocent_end_time] )
    fillTable( page, ['end_time_ns', inj.geocent_end_time_ns] )    
    fillTable( page, ['distance', '%.1f' % inj.distance] )
    fillTable( page, ['eff_dist_h','%.1f' %  inj.eff_dist_h] )
    fillTable( page, ['eff_dist_l','%.1f' %  inj.eff_dist_l] )
    fillTable( page, ['eff_dist_v','%.1f' %  inj.eff_dist_v] )
    fillTable( page, ['eff_dist_g','%.1f' %  inj.eff_dist_g] )  
    fillTable( page, ['playground','%s' %  InspiralUtils.isPlayground(inj)] )    
    page.add('</table></td>')
    
    # print infos to screen as well
    self.print_inj( inj,  injID)

    # now we can retrieve all the other files for this particular
    # missed injection
    foundDict = {}
    for stage, cache in self.triggerCache.iteritems():

      if self.exttrig:
        trigCache = lal.Cache(
          [c for c in cache if ('_'+str(injID)+'-' in os.path.basename(c.url)\
                                and inj.geocent_end_time in c.segment)])
      else:
        trigCache = lal.Cache(
          [c for c in cache if inj.geocent_end_time in c.segment])

      # check if the pfnlist is empty. `
      if trigCache.pfnlist()==[]:
        print >>sys.stderr, "Error: No files found for stage %s in the "\
              "cache for ID %s and time %d; probably mismatch of a "\
              "pattern in the options. " % \
              ( stage, injID, inj.geocent_end_time)
        continue

      # now create several plots
      if 'TMPLTBANK' in stage:
        horizon = self.expectedHorizonDistance( \
          trigCache.sieve(description=description).pfnlist(), \
          inj, ifo, self.number )
        foundDict[stage] = horizon
      elif 'INSPIRAL' in stage:
        found = self.investigateInspiral( \
          trigCache.sieve(description=description).pfnlist(), \
          inj, ifo, stage, self.number )
        foundDict[stage] = found
      elif 'THINCA' in stage:
        if 'FIRST' in stage:
          found = self.investigateThinca( \
            trigCache.sieve(description=description).pfnlist(), \
            inj, ifo, stage, self.number)     
        elif 'SECOND' in stage:
          # by default no "CAT_1" is associated to cat 1 times
          found = self.investigateThinca( \
            trigCache.sieve(description=description, exact_match=True).pfnlist(),\
            inj, ifo, stage, self.number )
        foundDict[stage]=found
      elif 'TRIGBANK' in stage:
        found = self.investigateTrigbank( \
          trigCache.sieve(description=description).pfnlist(), \
          inj, ifo, stage, self.number )
        foundDict[stage]=found
      else:
        print >>sys.stderr, "Error: Unknown pipeline stage ", stage
        
    ## print out the result for this particular injection
    page.add('<td><table border="2" >')
    fillTable( page, ['<b>step','<b>F/M'] )
    for stage in self.stageLabels:
      if foundDict.has_key( stage ):
        if "TMPLTBANK" in stage:
          for ifo in foundDict[stage].keys():
            text = "Expected SNR is " + ifo
            this_snr = "%.2f" % foundDict[stage][ifo]
            fillTable( page, [ text,  this_snr]) 
      
        elif foundDict[stage]:
          foundIFOs=''
          if "INSPIRAL" in stage:
            foundIFOs=' in '
            for i in foundDict[stage]:
              foundIFOs+=i+' '

          fillTable( page, [ stage,  'FOUND'+foundIFOs])          
        else:
          fillTable( page, [ stage,  '<font color="red">MISSED'])           
    page.add('</table>')
    page.add('</td></tr></table><br><br>')

    ## add the pictures to the webpage
    for stage, fname in self.fnameDict.iteritems():
      if stage!="TMPLTBANK":
        page.a(extra.img(src=[fname], width=400, \
                       alt=fname, border="2"), title=fname, href=[ fname ])
      
    # and write the html file
    htmlfilename = self.opts.prefix + "_"+ifo+"_followup_"+str(self.number) +\
                         self.opts.suffix+'.html'
    file = open(self.opts.output_path+htmlfilename,'w')      
    file.write(page(False))
    file.close()
    
    # supply the output
    return htmlfilename


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

