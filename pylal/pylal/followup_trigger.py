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

__Id__ = "$Id$"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

import os
import sys
import copy
from math import sqrt, pi
import subprocess
import tempfile

import matplotlib
matplotlib.use('Agg')
import pylab 
## try:
##   set
## except NameError:
##   from sets import Set as set

from pylal import SnglInspiralUtils
from pylal import InspiralUtils
from pylal import SimInspiralUtils
from pylal import CoincInspiralUtils
from pylal import SearchSummaryUtils
from pylal import grbsummary
from glue import lal
from glue import markup
from glue import segments
from glue import segmentsUtils
from glue.markup import oneliner as extra
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw.utils import ligolw_add

import numpy


##########################################################
class FollowupTrigger:
  """
  This defines a class for following up a trigger and to create
  the timeseries of found triggers in each stage of the pipeline.

  Usage:

  # first need to initialize the class
  followup = followup_trigger.FollowupTrigger(cache, opts)

  # later, one can do a followup of different objects, like an injection,
  #   a coinc. trigger or a single trigger:
  htmlname1 = followup.from_coinc(coinc)
  htmlname2 = followup.from_sngl(sngl_inspiral)
  htmlname3 = followup.from_missed(missed_inj)
  htmlname4 = followup.from_found(found_inj)

  # In each case the path to the created html file is returned. 
  # In the first call a CoincInspirals table is expected, a SngleInspiral
  # table in the second case and a SimInspiral table in the two last. 
  """


  # -----------------------------------------------------
  def __init__(self, cache, opts, use_injections = True):
    """
    Initialize this class and sets up all the cache files.
    @param cache: The cache of all files
    @param opts: The 'opts' structure from the main code
    @param use_injections: Specifying if injections are being used
                           (if no injections are used and this is set to
                           False, it speeds up the whole initalization...)
    """
    pylab.rcParams.update({'text.usetex': False})

    # Check all the required options
    option_list = ['verbose','followup_exttrig','output_path',\
                   'followup_time_window','prefix',\
                   'suffix','figure_resolution','user_tag']
    for option in option_list:
      if not hasattr(opts, option):
        raise "Error: The following parameter is required in the "\
              "opts structure for followup_trigger: ", option
    

    # setting the color definition and the stages of the pipeline
    self.colors = {'H1':'r','H2':'b','L1':'g','V1':'m','G1':'c'}
    self.stageLabels = ['INSPIRAL_FIRST', 'THINCA_FIRST',\
                        'INSPIRAL_SECOND', 'THINCA_SECOND']
    self.orderLabels = copy.deepcopy(self.stageLabels)
    self.orderLabels.extend( [ 'THINCA_SECOND_CAT_1','THINCA_SECOND_CAT_2', \
                               'THINCA_SECOND_CAT_3','THINCA_SECOND_CAT_4'] )

    # set arguments from the options
    self.opts = opts
    self.cache = cache
    self.tag = opts.user_tag
    self.verbose = opts.verbose
    self.exttrig = opts.followup_exttrig
    self.output_path = opts.output_path
    self.time_window = opts.followup_time_window
    self.sned = None
    
    # Value for the injection window. This value might be taken from 
    # a processParams table (see function get_injection_window)
    self.get_injection_window()

    # initialize a list of images created
    self.fname_list = []

    # counter for the followups
    self.number = 0

    # setting for each time to followup:
    self.followup_time = None
    self.injection_id = None  # only needed if time is an injection
    self.flag_followup = False

    if self.verbose:
      print "\nStarting initializing the Followup class..."
      
    # splitting up the cache for the different stages 
    self.trigger_cache = {}
    for stage in self.stageLabels:
      pattern = stage
      self.trigger_cache[stage] = self.cache.sieve(description=pattern)
      if self.opts.verbose:
        print "%d files found for stage %s" % (len(self.trigger_cache[stage]),\
                                               stage)


    # generating a dictionary for injection followups
    self.injections = dict()
    if use_injections:
      if self.verbose:
        print "Creating the injection cache..."

      self.injection_cache = self.cache.sieve(description = "INJECTION").\
                             sieve(ifos='HL')
      
      for file, entry in zip(self.injection_cache.pfnlist(), \
                             self.injection_cache):
        injection_id = self.get_injection_id(cache_entry = entry)
        self.injections[injection_id] = SimInspiralUtils.\
                                          ReadSimInspiralFromFiles( \
                                           [file], verbose=False )

    # read the veto files
    self.read_veto_files()
    
    if self.verbose:
      print "Initializing the Followup class done..."

  # -----------------------------------------------------
  def set_sned(self, executable):
    """
    Sets the sned flag to 'flag'
    @param flag: the executable used to calculate the sned
    """
    self.sned = executable

  # -----------------------------------------------------
  def execute_sned(self, inj):
    """
    Makes an external call to lalapps_sned to recalculate
    the effective distances
    @param inj: the injection to be cvonverted
    @return: the recalculated injection
    """

    
    file1 = 'test1.xml'
    file2 = 'test2.xml'
    
    # write out a dummy file
    grbsummary.write_rows([inj], lsctables.SimInspiralTable, file1)
  
    # call command for lalapps_sned
    command = self.sned+' --f-lower 40.0  --inj-file '+file1+\
              '  --output '+file2
    
    os.system(command)

    # read in the 'converted' injection
    doc = ligolw_add.ligolw_add(ligolw.Document(), [file2])
    inj_sned = lsctables.getTablesByType(doc, lsctables.SimInspiralTable)

    # return a single SimInspiral table
    return inj_sned[0][0]

  # -----------------------------------------------------
  def setTag(self, tag):
    """
    Setting a tag (non-conformal naming because of backwards compatibality)
    @param tag: the tag to be set
    """
    self.set_tag(tag)

  def set_tag(self, tag):
    """
    Setting a tag
    @param tag: well, its the tag!
    """
    self.tag = tag
    
      
  # -----------------------------------------------------
  def get_injection_window(self):
    """
    Extracting the length of the used injection_window from
    any 'FOUND' file in the cache.
    """

    # default value
    self.injection_window = 0.050

    # get the process params table from one of the COIRE files
    found_cache = self.cache.sieve(description = "FOUND")
    if len(found_cache)==0:
      # obviously no injections are being used. setting this window to zero
      self.injection_window = 0
      print "INFO: No FOUND files found, so setting the injection window"\
            " to zero."
      return
      
    coire_file = found_cache.checkfilesexist()[0].pfnlist()[0]
    try:
      doc = SearchSummaryUtils.ReadTablesFromFiles([coire_file],\
                                                   [lsctables.ProcessParamsTable])
      process_params = table.get_table(doc, lsctables.ProcessParamsTable.\
                                       tableName)
    except IOError: 	    
      sys.stderr.write("ERROR (IOError) while reading process_params table from"\
                       " file %s. Does this file exist and does it contain"\
                       " a search_summary table?\n" %(coire_file))
      raise 	 
    except AttributeError: 	 
      sys.stderr.write("ERROR (AttributeError:) while reading process_params"\
                       " table from file %s. Is the version of "\
                       "SearchSummaryUtils.py at least 1.5? Seems you have %s.\n" \
                       %(coire_file, SearchSummaryUtils.__version__))
      raise
    except:
      raise "Error while reading process_params table from file: ", coire_file

    # and retrieve the time window from this file
    found_flag = False
    for tab in process_params:
      if tab.param=='--injection-window':
        found_flag = True
        self.injection_window = float(tab.value)/1000.0
    if not found_flag: 	 
      sys.stderr.write("WARNING: No entry '--injection-window' found in file %s"\
                       "Value used is %.1f ms. If incorrect, please change file at %s\n" %\
                       (coire_file, 1000.0*self.injection_window, __file__))

    # set the parameter
    if self.verbose:
      print "Injection-window set to %.0f ms" % (1000*self.injection_window)

  
  # -----------------------------------------------------
  def get_injection_id(self, filename=None, url=None, cache_entry=None):
    """
    Extracting the injection ID from the filename, using
    the mechanism as used in lalapps_path2cache. You can specify the filename
    itself, the url or the cache entry. You must not specify more than one input!
    The injection-ID is calculated in the following way (exttrig only):

    The code expects the INSPIRAL and THINCA files in the following scheme (example):
      PREFIX-TAGPART_injections32_77-GPS-DURATION.xml
    The number of the injection run is extracted (INJRUN) as well as the following
    number (INJNUMBER). The injection ID is then calculated as:
       INJID = 100000*INJRUN + INJNUMBER
    so for this example the injectionj ID is 3200077. 
    
    @param filename: filename from which the injection ID is extracted
    @param url:  url from which the injection ID is extracted
    @param cache_entry: cache entry from which the injection ID is extracted
    """
    
    # Check that only one input is given
    if cache_entry:
      if filename and url:
        raise "Error in function get_injection_id: Only one input should be "\
              "specified. Now 'cache_entry' and another variable is specified. Check the code."
    
    if cache_entry is None:
      if filename and url:
        raise "Error in function get_injection_id: Only one input should be "\
              "specified. Now 'filename' and 'url' is specified. Check the code."
      
      if filename:
        path, filename = os.path.split(filename.strip())
        url = "file://localhost%s" % os.path.abspath(os.path.join(path, filename))

      try:
        cache_entry = lal.CacheEntry.from_T050017(url)
      except ValueError, e:
        raise "Error while extracting injection ID from file ", filename

    # split the expression into different sub-pieces
    pieces = cache_entry.description.split('_')
    if self.exttrig:

      # its easy for the exttrig case
      injection_id = pieces[-2]+'_'+pieces[-1]
    else:

      # but need to check for the appearance of the CAT suffix else
      index = 0
      for ind, piece in enumerate(pieces):
        if 'CAT' in piece:
          index = ind          
      injection_id = pieces[index-1] 
        
    return injection_id
  
  # -----------------------------------------------------  
  def find_injection_id(self, injection):
    """
    Find the injection-ID corresponding to this particular injection.
    @param injection: the injection object
    @return: the injection ID
    """

    if self.injections:
      # injection_id: the injection ID (a number or a string)
      # group_inj: a list of SimInspiral tables
      for injection_id, group_inj in self.injections.iteritems():
        for inj in group_inj:
          if injection.geocent_end_time==inj.geocent_end_time and \
                 injection.geocent_end_time_ns==inj.geocent_end_time_ns:
            return injection_id
          
      raise "No injection ID found for the above particular missed Injection "

    else:

      # cache not presearched, so searching for this particular injection
      if self.verbose:
        print "INFO: Searching for the injection at time %d.%d" % \
              (injection.geocent_end_time, injection.geocent_end_time_ns)
      injection_cache = self.cache.sieve(description = "INJECTION").\
                        sieve(ifos='HL')
      
      for file, entry in zip(injection_cache.pfnlist(), injection_cache):
        injection_id = self.get_injection_id(cache_entry = entry)
        sims = SimInspiralUtils.ReadSimInspiralFromFiles( [file], verbose=False )
        
        for sim in sims:
          # searching just by comparing the times...
          if injection.geocent_end_time==sim.geocent_end_time and \
                 injection.geocent_end_time_ns==sim.geocent_end_time_ns:
            if self.verbose:
              print "... found it: injection_id = ", injection_id
            return injection_id
      

    return None

  # -----------------------------------------------------
  def read_veto_files( self ):
    """
    Reads the veto segments given by veto-files (if any)
    """
    self.vetodict = dict()

    # loop over the IFO names
    for ifoName in self.colors.keys():

      self.vetodict[ifoName]=None
      # create the attribute name and check if it exists
      attributeName = 'followup_vetofile_'+ifoName.lower()
      if hasattr( self.opts, attributeName):

         # get the filename
         filename = getattr( self.opts, attributeName )
         if filename:
           self.vetodict[ifoName] = segmentsUtils.fromsegwizard(open(filename))

  # -----------------------------------------------------
  def reset( self ):
    """
    Resets the counting number for the time-series plots generated.
    """
    self.number=0

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
          ( inj.mass1, inj.mass2, inj.geocent_end_time, inj.geocent_end_time_ns,\
            inj.distance, inj.eff_dist_h, inj.eff_dist_l )

  # ----------------------------------------------------
  def save_plot( self, stage ):
    """
    Saves the plots and store them in a seperate fname_list.
    @param stage: the stage this plot belongs to (e.g. INSPIRAL, THINCA,...)
    """
    fname = 'Images/'+self.opts.prefix + "_"+self.tag+"_map-"+\
            stage+"-"+str(self.number) +self.opts.suffix+'.png'
    fname_thumb = InspiralUtils.\
                  savefig_pylal( filename = self.output_path+fname,\
                                 doThumb = True, 
                                 dpi_thumb = self.opts.figure_resolution)
    
    self.fname_list.append( fname ) 
    return fname
 

  # -----------------------------------------------------  
  def get_time_trigger( self, trig ):
    """
    This is a helper function to return a GPS time as one float number
    @param trig: a sngl_inspiral table entry
    """
    return float(trig.end_time) + float(trig.end_time_ns) * 1.0e-9

  # -----------------------------------------------------
  def get_sim_time(self, sim, ifo = None):
    """
    This is a helper function to return a GPS time as one float number
    for a certain IFO. If no IFO is specified the injected geocentric
    time is returned.
    @param sim: a sim_inspiral table entry
    @param ifo: the IFO for which we want the sim time
    """
    
    time=0
    nano=0

    if not ifo:
      time = sim.geocent_end_time
      nano = sim.geocent_end_time_ns
    if ifo:
      time = getattr(sim, ifo[0].lower()+'_end_time' )
      nano = getattr(sim, ifo[0].lower()+'_end_time_ns' )    

    return  float(time) + float(nano) * 1.0e-9


  # -----------------------------------------------------
  def is_veto(self, time_trigger, ifo):
    """
    This function checks if there is a veto at the time 'timeTrigger'
    for the IFO 'ifo'.
    @param time_trigger: The time to be investigated
    @param ifo: The name of the IFO to be investigated
    """
    
    return iterutils.any(time_trigger in seg for seg in self.vetodict[ifoName])

  # -----------------------------------------------------
  def put_text(self, text):
    """
    Puts some text into an otherwise empty plot.
    @param text: text to put in the empty plot
    """
    
    newText = ''
    for i in range( int(len(text)/60.0)+1):
      newText+=text[60*i:60*i+60]+'\n'
    pylab.figtext(0.15,0.15, newText)
 
  # -----------------------------------------------------
  def create_timeseries(self, trigger_files, stage, number):
    """
    Investigate inspiral triggers and create a time-series
    of the SNRs around the injected time
    @param trigger_files: List of files containing the inspiral triggers
    @param stage:        the name of the stage (FIRST, SECOND)
    @param number:       the consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose:
      print "Processing INSPIRAL triggers from files ", trigger_files
      
    sngls = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
              trigger_files , verbose=False)

    # create a figure and initialize some lists
    fig=pylab.figure()
    foundSet = set()
    loudest_details = {}
    no_triggers_found = True
    
    if sngls is None:
      self.put_text( 'No sngl_inspiral triggers in %s' % str(trigger_files))

    else:
      # create the small and large segments
      seg_small =  segments.segment(self.followup_time - self.injection_window, \
                                    self.followup_time + self.injection_window)
      seg_large =  segments.segment(self.followup_time - self.time_window, \
                                    self.followup_time + self.time_window)

      # create coincidences for THINCA stage
      coincs = None
      if 'THINCA' in stage:
        coincs = CoincInspiralUtils.coincInspiralTable(sngls, \
                   CoincInspiralUtils.coincStatistic("snr"))
        selected_coincs = coincs.vetoed(seg_small)
      
      # loop over the IFOs 
      for ifo in self.colors.keys():

        # get the singles for this ifo
        sngls_ifo = sngls.ifocut(ifo)

        # select a range of triggers
        selected_large = sngls_ifo.vetoed(seg_large)
        time_large = [ float(sel.get_end()) - self.followup_time \
                      for sel in selected_large ]

        selected_small = sngls_ifo.vetoed(seg_small)
        time_small = [ float(sel.get_end()) - self.followup_time \
                      for sel in selected_small ]

        # use the set of selected coincident triggers in the THINCA stages
        if coincs:
          selected_small = selected_coincs.cluster(2 * self.injection_window).getsngls(ifo)
          time_small = [ float(sel.get_end())-self.followup_time \
                        for sel in selected_small ]
          
        # skip if no triggers in the large time window
        if len(time_large)==0:
          continue
        no_triggers_found = False

        # add IFO to this set; the injection is found for this IFO-stage
        if len(time_small)>0:
          foundSet.add(ifo)                  

          # record details of the loudest trigger
          loudest = selected_small[selected_small.get_column('snr').argmax()]
          loudest_details[ifo] = {}
          loudest_details[ifo]["snr"] = loudest.snr
          loudest_details[ifo]["mchirp"] = loudest.mchirp
          loudest_details[ifo]["eff_dist"] = loudest.eff_distance
          loudest_details[ifo]["chisq"] = loudest.chisq
          loudest_details[ifo]["timeTrigger"] = float(loudest.get_end())

        # plot the triggers
        pylab.plot( time_large, selected_large.get_column('snr'),\
              self.colors[ifo]+'o', label="_nolegend_")
        pylab.plot( time_small, selected_small.get_column('snr'), \
              self.colors[ifo]+'s', label=ifo)

        # highlight any veto
        # FOR NOW: COMMENTED OUT...
        ## self.highlight_veto(self.followup_time, seg_large, ifo, ylims)  

      # draw the injection times and other stuff
      if no_triggers_found:
        self.put_text( 'No triggers/coincidences found within time window')
        
      ylims=pylab.axes().get_ylim()
      pylab.plot([0,0], ylims, 'g--', label="_nolegend_")
      pylab.plot([-self.injection_window, -self.injection_window], ylims, 'c:',\
            label="_nolegend_")
      pylab.plot([+self.injection_window, +self.injection_window], ylims, 'c:',\
            label="_nolegend_")

      # save the plot
      pylab.grid(True)
      pylab.legend()

    ylims=pylab.axes().get_ylim()
    pylab.axis([-self.time_window, +self.time_window, ylims[0], ylims[1]])
    pylab.xlabel('time [s]')
    pylab.ylabel('SNR')
    pylab.title(stage+'_'+str(self.number))    
    fname = self.save_plot( stage )
    pylab.close(fig)

    result = {'filename':fname, 'foundset':foundSet, 'loudest_details':loudest_details}
    return result


  # -----------------------------------------------------    
  def select_category(self, trigger_files, category):
    """
    Return a trigger list that contains only files for the choosen category.
    @param trigger_files : a list of trigger file names
    @param category: a category number
    @return: a sub list of filename corresponding to the category requested
    """

    # there are two different labels used to denote
    # the categories. THIS NEEDS TO BE UNIFIED
    cat1 = 'CAT_'+str(category)
    cat2 = 'CATEGORY_'+str(category)
    
    if category==1:
      # Category 1 files might not be labelled with a 'CAT_1' suffix.
      # So, for now, all files NOT containing the expression
      # 'CAT' in the filename are supposed to be CAT_1 files
      new_list = [file for file in trigger_files \
                  if 'CAT' not in file or cat1 in file or cat2 in file]
                     
    else:
      cat = 'CAT_'+str(category)
      new_list = [file for file in trigger_files if cat1 in file\
                  or cat2 in file]      
          
    return new_list

  # -----------------------------------------------------  
  def fill_table(self, page, contents):
    """ 
    Fills contents in a html table
    @param page: the pagfe object describing a html page
    @contents: the contents of the next table item
    """

    page.add('<tr>')
    for content in contents:
      page.add('<td>')
      page.add( str(content) )
      page.add('</td>')
    page.add('</tr>')

 
  # --------------------------------------------
  def create_table_inj(self, inj):
    """
    Creates the first table containing basic properties
    of the injection which is followed up.
    @param inj: an injection table
    """

    if self.sned:
      inj = execute_sned(inj)

    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup injection #"+str(self.number))
    page.add('<table border="2" >')          
    self.fill_table( page, ['<b>parameter','<b>value'] )
    self.fill_table( page, ['Number', self.number] )
    self.fill_table( page, ['inj ID', self.injection_id] )
    self.fill_table( page, ['mass1', '%.2f'% inj.mass1] )
    self.fill_table( page, ['mass2', '%.2f'% inj.mass2] )
    self.fill_table( page, ['mtotal', '%.2f' % (inj.mass1+inj.mass2)] )
    self.fill_table( page, ['mchirp', '%.2f' % (inj.mchirp)] )
    self.fill_table( page, ['end_time', inj.geocent_end_time] )
    self.fill_table( page, ['end_time_ns', inj.geocent_end_time_ns] )    
    self.fill_table( page, ['distance', '%.1f' % inj.distance] )
    self.fill_table( page, ['eff_dist_h','%.1f' %  inj.eff_dist_h] )
    self.fill_table( page, ['eff_dist_l','%.1f' %  inj.eff_dist_l] )
    self.fill_table( page, ['eff_dist_v','%.1f' %  inj.eff_dist_v] )
    self.fill_table( page, ['eff_dist_g','%.1f' %  inj.eff_dist_g] )  
    self.fill_table( page, ['playground','%s' %  InspiralUtils.isPlayground(inj)] )    
    page.add('</table></td>')
    page.hr()
    
    return page
  
  # --------------------------------------------
  def create_table_sngl(self, trig):
    """
    Creates the first table containing basic properties
    of the trigger which is followed up.
    @param trig: an sngl_inspiral table
    """
    
    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup trigger #"+str(self.number))
    page.add('<table border="2" >')          
    self.fill_table( page, ['<b>parameter','<b>value'] )
    self.fill_table( page, ['Number', self.number] )
    self.fill_table( page, ['inj ID', self.injection_id] )
    self.fill_table( page, ['mass1', '%.2f'% trig.mass1] )
    self.fill_table( page, ['mass2', '%.2f'% trig.mass2] )
    self.fill_table( page, ['mtotal', '%.2f' % (trig.mass1+trig.mass2)] )
    self.fill_table( page, ['mchirp', '%.2f' % (trig.mchirp)] )
    self.fill_table( page, ['end_time', trig.end_time] )
    self.fill_table( page, ['end_time_ns', trig.end_time_ns] )    
    self.fill_table( page, ['eff_distance', '%.1f' % trig.eff_distance] )
    page.add('</table></td><br>')
    page.hr()
    
    return page
  
  # --------------------------------------------
  def create_table_coinc(self, coinc):
    """
    Creates the first table containing basic properties
    of the coincidence which is followed up.
    @param coinc: an CoincInspiral table
    """
    
    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup trigger #"+str(self.number))
    page.add('<table border="2">')

    for ifo in ['H1','H2','L1','V1','G1']:
      if hasattr(coinc,ifo):
        trig = getattr(coinc,ifo)
        page.add('<td><table border="2" >')        
    
        self.fill_table( page, ['<b>parameter','<b>'+ifo] )
        self.fill_table( page, ['Number', self.number] )
        self.fill_table( page, ['inj ID', self.injection_id] )
        self.fill_table( page, ['SNR', trig.snr] )
        self.fill_table( page, ['ChiSq', trig.chisq] )
        self.fill_table( page, ['RSQ', trig.rsqveto_duration] )                        
        self.fill_table( page, ['Mass1', '%.2f'% trig.mass1] )
        self.fill_table( page, ['Mass2', '%.2f'% trig.mass2] )
        self.fill_table( page, ['Mtotal', '%.2f' % (trig.mass1+trig.mass2)] )
        self.fill_table( page, ['Mchirp', '%.2f' % (trig.mchirp)] )
        self.fill_table( page, ['end_time', trig.end_time] )
        self.fill_table( page, ['end_time_ns', trig.end_time_ns] )    
        self.fill_table( page, ['eff_distance', '%.1f' % trig.eff_distance] )
        page.add('</table></td>')                

    
    page.add('</table><br>')
    page.hr()
    
    return page
  
  # --------------------------------------------
  def create_table_time(self, trigger_time):
    """
    Creates the first table containing the time
    of the followup
    @param trigger_time: well, just the trigger time
    """
    
    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup time around GPS "+str(trigger_time) )
    page.add('<table border="2" >')          
    self.fill_table( page, ['Number', self.number] )
    self.fill_table( page, ['Time', trigger_time] )    
    page.add('</table></td><br>')
    page.hr()

    return page


  # --------------------------------------------
  def add_table_followup(self, page, invest_dict):
    """
    Adds the table containing specific information on the loudest
    trigger found in the followup region for each IFO.
    @param page: the html page to which to add the table
    @param invest_dict: dictionary containing the stage results
    """
    
    ## print out the result for this particular injection
    page.add('<td><table border="2" >')
    self.fill_table( page, ['<b>step','<b>F/M', '<b>Rec. SNR', '<b>Rec. mchirp', \
                      '<b>Rec. eff_dist', '<b>Rec. chisq', '<b>Veto ON/OFF'] )

    # loop over the stages and create the table with
    # the various data in it (when available)
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]

        # Fill in the details of the loudest found coinc.
        found_ifo=''
        loudest_snr=''
        loudest_mchirp=''
        loudest_eff_dist=''
        loudest_chisq=''
        veto_onoff=''

        # add all the IFO's for this coincident
        for ifo in result['foundset']:
          found_ifo += ifo+' '
          
          # Parameters of the loudest trigger, taken from the
          # 'loudest-details' dictionary, created in 'create_timeseries'
          loudest_snr += ifo + ': ' + str(result['loudest_details'][ifo]['snr'])+'<br>'
          loudest_mchirp += ifo + ': ' + str(result['loudest_details'][ifo]['mchirp'])+'<br>'
          loudest_eff_dist += ifo + ': ' + str(result['loudest_details'][ifo]['eff_dist'])+'<br>'
          loudest_chisq += ifo + ': ' + str(result['loudest_details'][ifo]['chisq'])+'<br>'
          
          # Check whether some of the ifo times is vetoed
          time_trigger = float(result['loudest_details'][ifo]['timeTrigger'])
          if self.vetodict[ifo]:
            veto = self.is_veto(time_trigger, ifo)
            veto_txt = 'OFF'
            if veto:
              veto_txt = 'ON'              
            veto_onoff+=ifo+': '+veto_txt+'<br>'
          else: 
            veto_onoff+=ifo+': No info<br>'

        # Fill the table whether something is found or not
        if len(result['foundset'])>0:
          self.fill_table( page, [ stage,  'FOUND in '+found_ifo, \
                                   'loudest<br>'+loudest_snr, \
                                   'loudest<br>'+loudest_mchirp, \
                                   'loudest<br>'+loudest_eff_dist,\
                                   'loudest<br>'+loudest_chisq, veto_onoff])
        else:
          self.fill_table( page, [ stage,  '<font color="red">MISSED'])
          
    page.add('</table>')
    page.add('</td></tr></table><br><br>')

    return page


  # -----------------------------------------------------  
  def from_coinc(self, coinc, ifo = None, more_infos = False, \
                 injection_id = None):
    """
    Creates a followup page from a coincident trigger.
    @param coinc: the coincidence to be followed up
    @param ifo: specifies the ifo to be used from the coinc.
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    if not ifo:
      print "WARNING (not bad): No IFO specified, using data from the "\
            "first IFO in the coincidence "
      for ifo_name in self.colors.keys():
        if hasattr(coinc, ifo_name):
          ifo = ifo_name
          break
    sngl = getattr(coinc, ifo)
    
    # set the time
    self.followup_time = float(sngl.get_end())
 
    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_coinc(coinc)
    self.flag_followup = more_infos

    # do the followup
    return self.followup(page)

  # -----------------------------------------------------
  def from_sngl(self, sngl, ifo = None, more_infos = False, injection_id = None):
    """
    Creates a followup page from a single trigger.
    @param sngl: the sngl trigger to be followed up
    @param ifo: NOT USED
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_sngl(sngl)
    self.flag_followup = more_infos

    # set the time
    self.followup_time = float(sngl.get_end())

    # do the followup
    return self.followup(page)
    
  # -----------------------------------------------------
  def from_missed(self, missed, ifo = None, more_infos = True, injection_id = None):
    """
    Creates a followup page from a missed injection.
    @param sngl: the missed injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    return self.from_injection(missed, ifo = ifo, more_infos = more_infos,\
                               injection_id = injection_id )    
    
  # -----------------------------------------------------
  def from_found(self, found, ifo = None, more_infos = False, injection_id = None):
    """
    Creates a followup page from a found injection.
    @param sngl: the found injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    return self.from_injection(found, ifo = ifo, more_infos = more_infos, \
                               injection_id = injection_id )
    
  # -----------------------------------------------------
  def from_injection(self, injection, ifo = None, more_infos = True, injection_id = None):
    """
    Creates a followup page from an injection.
    @param injection: the injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    # Set the injection ID if required
    if injection_id:
      self.injection_id = injection_id
    else:
      self.injection_id = self.find_injection_id(injection)

    # prepare the page
    page =  self.create_table_inj(injection)

    self.flag_followup = more_infos
    
    # set the time and do the followup
    self.followup_time = self.get_sim_time(injection, ifo)
    
    # do the followup
    return self.followup(page)

  
  # -----------------------------------------------------
  def from_time(self, trigger_time, ifo = None, more_infos = False, injection_id = None): 
    """
    Creates a followup page from a given time.
    @param trigger_time: the time to be followed up
    @param ifo: NOT USED
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    self.flag_followup = more_infos

    # prepare the page
    page =  self.create_table_time(trigger_time)

    # set the time
    self.followup_time = trigger_time
    self.injection_id = injection_id
    self.flag_followup = False    

    # do the followup
    return self.followup(page)
    

  # -----------------------------------------------------
  def followup(self, page):
    """
    Central followup procedure, finding corresponding files,
    generating the time-series and creating the output html files
    @param page: The head of the html page created with different informations
    @return: filename of the created html
    """
  
    # increase internal number:
    self.number+=1
    page.add('<br>')

    # loop over each stage
    invest_dict = {}    
    for stage, cache in self.trigger_cache.iteritems():

      # loop over each file in a stage
      trig_cache = lal.Cache()
      for c in cache:

        # check the time and the injection ID
        if self.followup_time in c.segment:
          if not self.injection_id or \
                 (self.injection_id and \
                  self.get_injection_id(url = c.url) == self.injection_id):
            trig_cache.append( c )
        
      # check if the pfnlist is empty. `
      file_list = trig_cache.pfnlist()
      if len(file_list)==0:
        print >>sys.stderr, "ERROR: No files found for stage %s in the "\
              "cache for ID %s and time %d; probably mismatch of a "\
              "pattern in the options. " % \
              ( stage, self.injection_id, self.followup_time)        
        continue

      # call the function to create the timeseries
      if 'THINCA_SECOND' in stage:
        # ... need to loop over the four categories
        for cat in [1,2,3,4]:          
          select_list=self.select_category(file_list, cat)
          if len(select_list)==0:
            print "WARNING (not that bad): No THINCA_SECOND files found for category ", cat
            continue          
          modstage = stage+'_CAT_' + str(cat)
          invest_dict[modstage] = self.create_timeseries(select_list,modstage, self.number)
      else:
        invest_dict[stage]=self.create_timeseries(file_list, stage, self.number)


    ## add some more followup if required
    if self.flag_followup:
      self.add_table_followup(page, invest_dict)

    ## add the pictures to the webpage
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]
      
        fname = result['filename']
        page.a(extra.img(src=[fname], width=400, \
                         alt=fname, border="2"), title=fname, href=[fname])


    ## add the version of this code
    page.add("<hr>")
    page.add("Figure(s) and data produced with " + __name__ + ", version " \
              + __version__)
        
    # and write the html file
    htmlfilename = self.opts.prefix + "_followup_"+str(self.number) +\
                         self.opts.suffix+'.html'
    file = open(self.opts.output_path+htmlfilename,'w')      
    file.write(page(False))
    file.close()

    # store html file in fname_list and return filename
    self.fname_list.append(htmlfilename)
    return htmlfilename

