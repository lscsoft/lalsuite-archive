#!/usr/bin/python

import sys, os, re, glob, exceptions
import numpy

from optparse import *

from glue import lal
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import segments
from glue import segmentsUtils


from pylal import SnglInspiralUtils
from pylal import CoincInspiralUtils
from pylal import InspiralUtils
from pylal import SimInspiralUtils

 
usage =  """Usage: %prog [options]

hanfordCheck \
--cache-file     CACHEFILE \
--gps-start-time GPSSTARTTIME \
--gps-end-time   GPSENDTIME \
--distance-cut   DISTANCECUT \
--h2-threshold   H2THRESHOLD \
--other-ifos     OTHERIFOS \
--user-tag       USERTAG \
--veto-file      VETOFILE \
--segs-file      SEGSFILE \
--skip-zero-lag \
--check-slides \
--verbose"""


#
# =============================================================================
#
#                     Parse Command Line
#
# =============================================================================
#

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage)
  parser.add_option("-c", "--cache-file", help="LAL cache of relevant files")
  parser.add_option("-S", "--gps-start-time", default="0", \
                           help="gps start time")
  parser.add_option("-E", "--gps-end-time", default="0", \
                           help="gps end time")
  parser.add_option("-V", "--veto-file", action="store", default=None, \
                           help="H1 SEGS file of relevant vetoes")
  parser.add_option("-f", "--segs-file", action="store", default=None, \
                           help="segs of time you want to check e.g H1H2L1")
  parser.add_option("-u", "--user-tag", help="User tag");
  parser.add_option("-H", "--h2-threshold", type="float", default=5.5, \
                     help="H2 snr threshold")
  parser.add_option("-o","--other-ifos", default="L1",\
                     help="e.g. L1, L1V1, V1 etc.")
  parser.add_option("-d", "--distance-cut", type="float", default=0.6,\
                     help="distance cut")
  parser.add_option("-e", "--enable-output", action="store_true", \
                     default=False, help="print html output" )
  parser.add_option("-z", "--skip-zero-lag", action="store_true", \
                     default=False, help="check zero lag triggers" )
  parser.add_option("-s", "--check-slides", action="store_true", \
                     default=False, help="check-slides" )
  parser.add_option("-v", "--verbose",action="store_true",\
                     default=False,help="print all information" )
  command_line = sys.argv[1:]
  (options,args) = parser.parse_args()

  return options, sys.argv[1:]


# ============================================================================
# -- get command line arguments
opts, args = parse_command_line()

# Must have one or the other!
if opts.skip_zero_lag:
  opts.check_slides = True

if opts.gps_start_time == 0 or opts.gps_end_time == 0:
  print "You must specify --gps-start/end-time."
  sys.exit()
times = segments.segment( opts.gps_start_time, opts.gps_end_time )

if not opts.veto_file:
  print "You must specify --veto-file."
  sys.exit()
vetoes = segmentsUtils.fromsegwizard( open( opts.veto_file ) )

if opts.segs_file:
  segs = segmentsUtils.fromsegwizard( open( opts.segs_file ) )


title = "hanford_check_" + opts.user_tag



##############################################################################
def get_coire_triggers( opts, ifo, pattern, times ):
  """
  Returns sngl triggers from coire files
  """

  trigcache = cache.sieve( ifos = ifo, description = pattern, segment = times )
  trigFiles = trigcache.checkfilesexist()[0].pfnlist()

  if len( trigFiles ) == 0:
    if opts.verbose:
      print "ERROR! Cannot find specified coire files. Check cache file", \
            "and/or specified ifos/injections/gps times.\n"
    sys.exit()

  if opts.verbose:
    print "Reading coire files..."
    if opts.verbose:
      for file in trigFiles:
        print file

  inspTriggers  = SnglInspiralUtils.ReadSnglInspiralFromFiles( trigFiles )
  coincTriggers = CoincInspiralUtils.coincInspiralTable( inspTriggers, \
                             CoincInspiralUtils.coincStatistic( "snr") )

  triggers = coincTriggers.getsngls( ifo )

  return triggers



##############################################################################
def check_triggers( opts, pattern, title ):

  html_table=[]
  html_move=[]

  # For trigger info
  doubles=[]
  triples=[]

  ifo = "H1"
  h1 = get_coire_triggers( opts, ifo, pattern, times)
  h1_events = list( h1.get_column( "event_id" ) )

  ifo = "H2"
  h2 = get_coire_triggers( opts, ifo, pattern, times)
  h2_events = list( h2.get_column( "event_id" ) )
  h2_times = list( h2.get_column( "end_time" ) )
  h2_snr = list( h2.get_column( "snr" ) )


  h2_index=[]
  for event in h2_events:
    if event not in h1_events:
      h2_index.append( h2_events.index( event ) )


  write_table_head( opts, html_table )

  for i in h2_index:

    columns = []

    # H2 Time and SNR
    columns.append( h2_times[i] )
    columns.append( h2_snr[i] )
  
    # H2 Range
    for j in xrange( len( h2SegStart ) ):
      if h2_times[i] >= h2SegStart[j] and h2_times[i] < h2SegEnd[j]:
        H2range = h2SegRange[j]  
    H1range = 0
    for j in xrange( len( h1SegStart ) ):
      if h2_times[i] >= h1SegStart[j] and h2_times[i] < h1SegEnd[j]:
        H1range = h1SegRange[j]  

    # Distance Cut
    H2dist = H2range * opts.h2_threshold / h2_snr[i]

    columns.append( H2dist )  
    columns.append( H2range )  
    columns.append( H1range )  

    cut = 2 * abs( H1range - H2dist ) / ( H1range + H2dist )
    columns.append( cut )

    # H1 Veto 
    if h2_times[i] in vetoes:
      columns.append( 'YES' )
      doubles.append( 1 )
    else:
      columns.append( 'NO' )
      triples.append( 1 )

    write_table_row( opts, columns, html_table )

  line = '</table>\n'
  html_table.append( line )



  # ===========================================
  # Trigger Info

  if opts.verbose:
    print "Number of H2" + opts.other_ifos + " trigers in H1H2"\
          + opts.other_ifos + " time = " + str( len( h2_index ) )
  if 'SLIDE' in pattern:
    line = '<h2> Trigger Info (Slides) </h2> \n'
  else:
    line = '<h2> Trigger Info (Zero Lag) </h2> \n'

  html_move.append( line )
  line = '<table cellpadding=5 border=1> \n'
  html_move.append( line )
  line = '<tr>\n'
  html_move.append( line )
  line = ' <td> Total H2' + opts.other_ifos + 'Triggers </td>\n'  
  html_move.append( line )
  line = ' <td>' + str( len ( h2_index ) ) + ' </td>\n'
  html_move.append( line )
  line = '</tr>\n'
  html_move.append( line )

  line = '<tr>\n'
  html_move.append( line )
  line = ' <td> Move to H2' + opts.other_ifos + ' time </td>\n'
  html_move.append( line )
  line = ' <td>' + str( len ( doubles ) ) + ' </td>\n'
  html_move.append( line )
  line = '</tr>\n'
  html_move.append( line )
  line = '<tr>\n'
  html_move.append( line )
  line = ' <td> Genuine H2L1 in H1H2' + opts.other_ifos + ' Time </td>\n'
  html_move.append( line )
  line = ' <td>' + str( len ( triples ) ) + ' </td>\n'
  html_move.append( line )
  line = '</tr>\n'
  html_move.append( line )
  line += '</table>\n'
  html_move.append( line ) 

  return html_table, html_move



# =============================================================================
#
#                    Generic html functions
#
# =============================================================================
#

def write_html_head( opts,  title ):
  """
  Creates a new html document wih the contents for the output
  """

  if opts.verbose:
    print "Initiating output web page..."

  # Begin html
  line = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
  line = line + '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> \n'
  html_head.append( line )
  line = '<html xmlns="http://www.w3.org/1999/xhtml"> \n'
  html_head.append( line )
  line = '<head> \n'
  html_head.append( line )
  line = '<meta http-equiv="Content-Type" content="text/html;'
  line = line + ' charset=UTF-8" /> \n'
  html_head.append( line )
  line = '<link media="all" href="style.css" type="text/css"'
  line = line + 'rel="stylesheet" /> \n'
  html_head.append( line )
  line = '<title>' + title + '</title> \n'
  html_head.append( line )
  line = '</head> \n'
  html_head.append( line )
  line = '<body> \n\n'
  html_head.append( line )
  line = '<h1>' + title + '</h1> \n' 
  html_head.append( line )
  line = '<h2> Arguments </h2> \n' 
  html_head.append( line )
  line = '<p> '
  html_head.append( line )
  line = '  <code> '
  html_head.append( line )
  for i in xrange( len( sys.argv ) ):
    line = line + sys.argv[i] + ' '
  line = line + ' \n'
  html_head.append( line )
  line = '  </code> \n'
  html_head.append( line )
  line = '</p> \n'
  html_head.append( line )



##############################################################################
def write_table_head( opts, html_table ):

  line = '<table cellpadding=5 border=1> \n'
  html_table.append( line )
  line = '<tr> \n'
  html_table.append( line )
  line = '  <th> H2 end time </th> \n'
  html_table.append( line )
  line = '  <th> H2 snr </th> \n'
  html_table.append( line )
  line = '  <th> H2 Distance </th> \n'
  html_table.append( line )
  line = '  <th> H2 Range </th> \n'
  html_table.append( line )
  line = '  <th> H1 Range </th> \n'
  html_table.append( line )
  line = '  <th> Distance Cut </th> \n'
  html_table.append( line )
  line = '  <th> H1 Vetoed ? </th> \n'
  html_table.append( line )
  line = '</tr> \n'
  html_table.append( line )

  return html_table



##############################################################################
def write_table_row( opts, columns, html_table ):
  """
  Assumes columns matches those set out in the write_table_head() function.
  """


  line = '<tr> \n'
  html_table.append( line )
  for i in xrange( len( columns ) ):
     
    # GPS time
    if i == 0:
      line = '  <td>' 
      line += str( columns[i] )
      line += '</td> \n'
    # SNR
    elif i in [1,2,3,4]:
      line = '  <td>' 
      line += '%.2f' % ( columns[i] )
      line += '</td> \n'
    # Distance Cut
    elif i == 5:
      if columns[i]  < opts.distance_cut:
        colour = "#00ff00"
      elif columns[i] > opts.distance_cut:
        colour = "#ff0000"
      else:
        colour = "#50ebec"
      line = '  <td bgcolor="' + colour +'">'
      line += '%.3f' % ( columns[i] )
      line += '</td>\n'
    # Veto?
    elif i == 6:
      if columns[i] == 'NO':
        colour = "#00ff00"
      elif columns[i] == 'YES':
        colour = "#ff0000"
      line = '  <td align="center" bgcolor="' + colour +'">'
      line += columns[i]
      line += '</td>\n'
    html_table.append( line )

  line = ' </tr> \n'
  html_table.append( line )

  return html_table



# =============================================================================
#
#                         Main                   
#
# =============================================================================
#

# ===========================================
# Declare arrays for output

html_head=[]
html_rhomax=[]
html_movers=[]
html_tables=[]
html_close=[]

write_html_head( opts, title )



# ===========================================
# Open cache file

cache = lal.Cache.fromfile( open( opts.cache_file ) )



# ===========================================
# Check ranges over entire time

pattern      = 'INSPIRAL_FIRST_FULL_DATA'
trigcache    = cache.sieve( "H2", description = pattern, segment = times )
trigFiles    = trigcache.checkfilesexist()[0].pfnlist()
h2Summ, null = InspiralUtils.readHorizonDistanceFromSummValueTable( trigFiles, opts.verbose )

h2SegStart =  h2Summ['H2'].getColumnByName('start_time').asarray()
h2SegEnd   =  h2Summ['H2'].getColumnByName('end_time').asarray()
h2SegRange =  h2Summ['H2'].getColumnByName('value').asarray()

trigcache    = cache.sieve( "H1", description = pattern, segment = times )
trigFiles    = trigcache.checkfilesexist()[0].pfnlist()
h1Summ, null = InspiralUtils.readHorizonDistanceFromSummValueTable( trigFiles, opts.verbose )

h1SegStart =  h1Summ['H1'].getColumnByName('start_time').asarray()
h1SegEnd   =  h1Summ['H1'].getColumnByName('end_time').asarray()
h1SegRange =  h1Summ['H1'].getColumnByName('value').asarray()

rhomax = [[],[]]

# Check over entire H2 time!
for i in xrange( len( h2SegStart ) ):
  for j in xrange( len( h1SegStart ) ):
    time = 0
    if h1SegStart[j] >= h2SegStart[i] and h1SegStart[j] < h2SegEnd[i]:
      if h1SegEnd[j] < h2SegEnd[i]:
        time = h1SegEnd[j] - h1SegStart[j]
      else:
        time = h2SegEnd[i] - h1SegStart[j]
    elif h1SegEnd[j] > h2SegStart[i] and h1SegEnd[j] <= h2SegEnd[i]:
      time = h1SegEnd[j] - h2SegStart[i]
    if time != 0:
      max = opts.h2_threshold * h2SegRange[i] / h1SegRange[j] 
      max*= ( 2 + opts.distance_cut ) / ( 2 - opts.distance_cut )
      rhomax[0].append( max )
      rhomax[1].append( time )

totaltime = sum( rhomax[1] )

line = '<h2> Range Information - not including vetoes </h2> \n'
html_rhomax.append( line )
line = '<h4> All H1H2 times </h4> \n'
html_rhomax.append( line )
line = '<ul><li> H2 SNR threshold = ' + str( opts.h2_threshold ) + '</li> \n'
html_rhomax.append( line )
line = '<li> Time analysed = ' + str( totaltime ) + ' seconds </li></ul> \n'
html_rhomax.append( line )
line = '<table cellpadding=5 border=1> \n'
html_rhomax.append( line )
line = '<tr>\n'
html_rhomax.append( line )
line = ' <th> &rho; Max </th>\n'
html_rhomax.append( line )
line = ' <th> time / s </th>\n'
html_rhomax.append( line )
line = ' <th> Percentage </th>\n'
html_rhomax.append( line )
line = '</tr>\n'
html_rhomax.append( line )

# Check threshold(s)
kmax = 5
for k in xrange( 0, kmax + 1 ):
  above = 0
  for i in xrange( len( rhomax[0]) ):
    if rhomax[0][i] > opts.h2_threshold + k * 0.1:
      above += rhomax[1][i]
  line = '<tr>\n'
  html_rhomax.append( line )
  line = ' <td> > ' + str( opts.h2_threshold + k * 0.1 ) + ' </td>\n'
  html_rhomax.append( line )
  line = ' <td> ' + str( above ) + ' </td>\n'
  html_rhomax.append( line )
  line = ' <td> ' + '%.2f'  % ( 100.0 * above / totaltime ) + ' % </td>\n'
  html_rhomax.append( line )
  line = '</tr>\n'
  html_rhomax.append( line )
line = '</table>\n'
html_rhomax.append( line )



# Check over selected segs!
if opts.segs_file:
  rhomax = [[],[],[],[]]

  for k in xrange( len( segs ) ):
    for i in xrange( len( h2SegStart ) ):
      theseTimes = 0
      if h2SegStart[i] >= segs[k][0] and h2SegStart[i] < segs[k][1]:
        if h2SegEnd[i] < segs[k][1]:
          theseTimes = [ h2SegStart[i], h2SegEnd[i] ] 
        else:
          theseTimes = [ h2SegStart[i], segs[k][1] ]
      elif h2SegEnd[i] > segs[k][0] and h2SegEnd[i] <= segs[k][1]:
        theseTimes = [ segs[k][0], h2SegEnd[i] ]
      
      if theseTimes != 0:	
        for j in xrange( len( h1SegStart ) ):
          time = 0
          if h1SegStart[j] >= theseTimes[0] and h1SegStart[j] < theseTimes[1]:
            if h1SegEnd[j] < theseTimes[1]:
              time = h1SegEnd[j] - h1SegStart[j]
              start =  h1SegStart[j]
              end =  h1SegEnd[j]
            else:
              time = theseTimes[1] - h1SegStart[j]
              start = h1SegStart[j] 
              end = theseTimes[1] 
          elif h1SegEnd[j] > theseTimes[0] and h1SegEnd[j] <= theseTimes[1]:
            time = h1SegEnd[j] - theseTimes[0]
            start = theseTimes[0]
            end = h1SegEnd[j]
          if time != 0:
            max = opts.h2_threshold * h2SegRange[i] / h1SegRange[j] 
            max*= ( 2 + opts.distance_cut ) / ( 2 - opts.distance_cut )
            rhomax[0].append( max )
            rhomax[1].append( time )
            rhomax[2].append( start )
            rhomax[3].append( end )

  totaltime = sum( rhomax[1] )

  line = '<h4> H1H2' + opts.other_ifos + ' Times: "<code> ' + opts.segs_file 
  line+= '</code>"</h4> \n'
  html_rhomax.append( line )
  line = '<ul><li> H2 SNR threshold = ' + str( opts.h2_threshold ) + '</li> \n'
  html_rhomax.append( line )
  line = '<li> Time analysed = ' + str( totaltime ) + ' seconds </li></ul> \n'
  html_rhomax.append( line )
  line = '<table cellpadding=5 border=1> \n'
  html_rhomax.append( line )
  line = '<tr>\n'
  html_rhomax.append( line )
  line = ' <th> &rho; Max </th>\n'
  html_rhomax.append( line )
  line = ' <th> time / s </th>\n'
  html_rhomax.append( line )
  line = ' <th> Percentage </th>\n'
  html_rhomax.append( line )
  line = '</tr>\n'
  html_rhomax.append( line )

  # Check threshold(s)
  kmax = 5
  for k in xrange( 0, kmax + 1 ):
    above = 0
    for i in xrange( len( rhomax[0]) ):
      if rhomax[0][i] > opts.h2_threshold + k * 0.1:
        above += rhomax[1][i]
        
    line = '<tr>\n'
    html_rhomax.append( line )
    line = ' <td> > ' + str( opts.h2_threshold + k * 0.1 ) + ' </td>\n'
    html_rhomax.append( line )
    line = ' <td> ' + str( above ) + ' </td>\n'
    html_rhomax.append( line )
    line = ' <td> ' + '%.2f'  % ( 100.0 * above / totaltime ) + ' % </td>\n'
    html_rhomax.append( line )
    line = '</tr>\n'
    html_rhomax.append( line )
  line = '</table>\n'
  html_rhomax.append( line )

  line = '<div style="position:relative; left:325px; top:-225px">\n'
  html_rhomax.append( line )
  line = '<h4> Corresponding Times (Where &rho; > ' + str( opts.h2_threshold ) 
  line+= ' )  </h4>\n'
  html_rhomax.append( line )
  line = '<ul>\n'
  html_rhomax.append( line )

  for i in xrange( len( rhomax[0] ) ):
    if rhomax[0][i] > opts.h2_threshold:
      line = '<li><code>' + str( rhomax[2][i] ) + ' - ' 
      line += str( rhomax[3][i] ) + ' = ' 
      line += str( rhomax[1][i] ) +' </code></li> \n'
      html_rhomax.append( line )

  line = '</ul>\n'
  html_rhomax.append( line )

  line = '</div>\n'
  html_rhomax.append( line )


# ===========================================
# Check zero lag and/or slide triggers

if not opts.skip_zero_lag:
  pattern = 'COIRE_SECOND_H1H2' + opts.other_ifos + '_FULL_DATA_CAT_3'
  table, move = check_triggers( opts, pattern, title )
  html_tables.append( table )
  html_movers.append( move )

if opts.check_slides:
  pattern = 'COIRE_SLIDE_SECOND_H1H2' + opts.other_ifos + '_FULL_DATA_CAT_3'
  table, move = check_triggers( opts, pattern, title )
  html_tables.append( table )
  html_movers.append( move )


# ===========================================
# Finish page

html_break = '<br/><br/><br/><hr/>\n'

line = '</body>\n'
html_close.append( line )
line = '</html>\n'
html_close.append( line )


if opts.enable_output:
  fout = open( title + ".html", "w" )
  for line in html_head:
    fout.write( line )
  fout.write( html_break )
  for line in html_rhomax:
    fout.write( line )
  fout.write( html_break )

  line = '<table border=0> \n'
  fout.write( line )
  line = '<tr> \n'
  fout.write( line )
  for i in xrange( len( html_movers ) ):
    line = '  <td valign="top"> '
    fout.write( line )
    for line in html_movers[i]:
      fout.write( line )
    for line in html_tables[i]:
      fout.write( line )
    line = '  </td> '
    fout.write( line )
  line = '</tr> \n'
  fout.write( line )
  line = '</table> \n'
  fout.write( line )

  for line in html_close:
    fout.write( line )
  fout.close()
