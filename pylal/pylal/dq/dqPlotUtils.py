#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import math,re
from numpy import arange

from glue.ligolw import table
from glue.segments import segment, segmentlist

from datetime import datetime
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal import date
from pylal.dq.dqTriggerUtils import def_get_time

import matplotlib
matplotlib.use( 'Agg' )
import pylab

from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides plotting routines for use in data quality investigations. All routines are written to work in as general a way as possible with ligolw tables and lsctables compatible columns. 
"""

# =============================================================================
# Function to plot before/after cumulative SNR histograms
# =============================================================================

def plot_trigger_hist( triggers, outfile, column='snr', segments=None,\
                       flag='unknowns', etg='Unknown',\
                       livetime=None, fill=False, logx=True,
                       cumulative=True, rate=True ):

  """
    Plot a histogram of the value in any column of the ligolw table triggers.
    If a glue.segments.segmentlist segments is given, the histogram is presented
    before and after removal of triggers falling inside any segment in the list.

    Arguments:

      triggers : [ SnglBurstTable | SnglInspiralTable | SnglRingdownTable ] 
        ligolw table containing triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      column : string
        valid column of triggers table to plot as histrogram
      segments : glue.segments.segmentlist
        list of segments with which to veto triggers
      flag : string
        display name of segmentlist, normally the name of the DQ flag
      etg : string
        display name of trigger generator, defaults based on triggers tableName
      livetime : [ float | int | LIGOTimeGPS ]
        span of time from which triggers and segments are valid, used to
        display histogram counts in terms of rate (Hz) for easy comparisons
      fill : [ True | False ]
        boolean option to fill below the histogram curves, default colors:
        red (vetoed), green (not vetoed).
      logx : [ True | False ]
        boolean option to display column (x) axis in log scale.
  """

  get_time = def_get_time( triggers.tableName )

  # calculate livetime
  if not livetime:
    times = [ get_time( t ) for t in triggers ] 
    start = min( times )
    end   = max( times )
    livetime = end-start 

  # generate vetoed trigger list: inspiral trigs have dedicated function 
  if not segments:
    segments = segmentlist()

  aftertriggers = table.new_from_template( triggers )
  aftertriggers.extend( [ t for t in triggers if get_time( t ) not in segments ] )

  # set up histogram data
  try:
    if logx:
      start_data = [math.log10( trig.__getattribute__( column ) )\
                    for trig in triggers]
      end_data   = [math.log10( trig.__getattribute__( column ) )\
                    for trig in aftertriggers]
    else:
      start_data = [trig.__getattribute__( column ) for trig in triggers]
      end_data   = [trig.__getattribute__( column ) for trig in aftertriggers]
  except KeyError:
    err = 'Column %s not found in %s.' % ( column,triggers.tableName )
    raise KeyError, err

  # generate histogram
  num_bins = 1000
  if len( triggers )>=1:

    if cumulative:
      cumulative = -1

    start_n,start_bins,start_p = pylab.hist( start_data, bins=num_bins,\
                                             range=( min( start_data ),\
                                                     max( start_data ) ),\
                                             histtype='stepfilled',\
                                             cumulative=cumulative,\
                                             facecolor='red',\
                                             visible=False )
    end_n,end_bins,end_p = pylab.hist( end_data, bins=num_bins,\
                                       range=( min( start_data ),\
                                               max( start_data ) ),\
                                       histtype='stepfilled',\
                                       cumulative=cumulative,\
                                       facecolor='green',\
                                       visible=False )

    # recalculate centre of bins ( in logscal, if required )
    bins = []
    for i in range( len( start_bins )-1 ):
      if logx:
        bins.append( 0.5 * ( math.pow( 10, start_bins[i] )+\
                             math.pow( 10, start_bins[i+1] ) ) )
      else:
        bins.append( 0.5 * ( start_bins[i]+start_bins[i+1] ) )

    # reset zero values to base ( so logscale doesn't break )
    if rate:
      livetime = float( livetime )
      base = 0.5/livetime
    else:
      base = 0.5
    def zero_replace( num, base ):
      if num==0:
        return base
      else:
        return num
    start_n = map( lambda n: zero_replace( n, base ), start_n )
    end_n   = map( lambda n: zero_replace( n, base ), end_n )

    # convert number to rate
    if rate:
      start_n = [n/livetime for n in start_n]
      end_n   = [n/livetime for n in end_n]

  else:
    bins = []
    start_n = []
    end_n = []

  # fix name for latex
  flag = flag.replace( '_','\_' )

  # customise plot appearance
  pylab.rcParams.update( {"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16,
                         "axes.titlesize": 20,
                         "axes.labelsize": 18,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20} )

  # plot data
  fig = pylab.figure( figsize=[12,6] )
  ax  = fig.gca()
  ax.loglog()

  if not fill:
    ax.plot( bins,start_n,'r',linewidth=2,label='Before vetoes' )
    if segments:
      ax.plot( bins,end_n,'g',linewidth=2,label='After vetoes' )
  if fill:
    ax.plot( bins,start_n,'r',linewidth=0,label='Before vetoes' )
    ax.fill_between( bins,base,start_n,color='r',edgecolor='k',linewidth=0.5 )
    if segments:
      ax.plot( bins,end_n,'g',linewidth=0,label='After vetoes' )
      ax.fill_between( bins,base,end_n,color='g',edgecolor='k',linewidth=0.5,\
                      alpha=0.9 )

  # figure sundries
  if segments:
    leg = ax.legend( loc='best' )
    for l in leg.get_lines():
      l.set_linewidth( 4 )
  if bins:
    ax.set_xlim( min( bins ),max( bins ) )
    ax.set_ylim( base,max( start_n )*1.01 )
  ax.set_xlabel( column.replace( '_','\_' ).title() )
  if rate and cumulative:
    ax.set_ylabel( 'Cumulative rate ( Hz )' )
  elif rate:
    ax.set_ylabel( 'Rate ( Hz )' )
  elif not rate and cumulative:
    ax.set_ylabel( 'Cumulative number' )
  elif not rate and not cumulative:
    ax.set_ylabel( 'Number' )
  tit = '%s triggers' % ( etg.replace( '_','\_' ) )
  if segments:
    tit += ' and %s segments' % ( flag )
  ax.set_title( tit )
  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )
  ax.set_axisbelow( True )
  fig.savefig( outfile )

# =============================================================================
# Function to plot SNR vs. time with veto segments
# =============================================================================

def plot_triggers( triggers, outfile, etg='Unknown',\
                   start=None, end=None, zero=None,\
                   segments=None, flag='unknown',\
                   xcolumn='time', ycolumn='snr', zcolumn=None,\
                   logx=False, logy=True, logz=True, ylim=None ):

  """
    Plots ycolumn against xcolumn for columns in given
    Sngl{Burst,Inspiral}Table object triggers, coloured by the zcolumn
    highlighting those entries falling inside one of the entries in the
    glue.segments.segmentlist object segments, if given. 

    'time' given as a column name is a special case, since s and ns times are
    stored separately in the SnglTable structures. In this case the
    trigger.get_xxxx() function is called.

    Arguments arguments:
    
      start : [ float | int | LIGOTimeGPS]
        GPS start time
      end : [ float | int | LIGOTimeGPS]
        GPS end time
      triggers : [ SnglBurstTable | SnglInspiralTable | SnglRingdownTable ] 
        ligolw table containing triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      etg : string
        display name of trigger generator, defaults based on triggers tableName
      start : [ float | int | LIGOTimeGPS]
        GPS start time
      end : [ float | int | LIGOTimeGPS]
        GPS end time
        list of segments with which to veto triggers
      flag : string
        display name of segmentlist, normally the name of the DQ flag
      xcolumn : string
        valid column of triggers table to plot on x-axis
      ycolumn : string
        valid column of triggers table to plot on y-axis
      zcolumn : string
        valid column of triggers table to use for colorbar. If not given, no
        colorbar will be used and all triggers will be blue.
      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
      logz : [ True | False ]
        boolean option to display z-axis (colorbar) in log scale.
  """

  get_time = def_get_time( triggers.tableName )

  # calculate livetime
  if not start or end:
    times = [ get_time( t ) for t in triggers ]
  if not start:
    start = int( math.floor( min( times ) ) )
  if not end:
    end   = int( math.ceil( max( times ) ) )

  if not zero:
    zero = start

  # sort triggers by z param
  if zcolumn:
    triggers.sort( key=lambda trig: trig.__getattribute__( zcolumn ),reverse=False )

  # apply veto segments if required
  if segments is not None:
    segments = segmentlist( segments )
  
    # set up vetoed/nonvetoed trigger lists, inspiral triggers have dicated func
    trigs =  table.new_from_template( triggers )
    vetotrigs = table.new_from_template( triggers )

    trigs.extend( t for t in triggers if get_time( t ) not in segments )
    vetotrigs.extend( t for t in triggers if get_time( t ) in segments )

  # or copy lists with empty vetoed triggers list
  else:
    trigs     = triggers
    vetotrigs = table.new_from_template( triggers )

  # set plot time unit whether it's used or not
  if ( end-start ) < 20000:
    t_unit = 60
  elif ( end-start ) >= 20000 and ( end-start ) < 604800:
    t_unit = 3600
  else:
    t_unit = 86400
  t_string = {60:'minutes',3600:'hours',86400:'days'}

  # set up plot lists
  notvetoed = {}
  vetoed = {}

  for col in [xcolumn,ycolumn]:
    try:
      # treat 'time as special case'
      if col=='time':
        notvetoed[col] = [ float( get_time(t) - zero ) / t_unit\
                           for t in trigs ]
        vetoed[col] = [ float( get_time(t) - zero ) / t_unit\
                        for t in vetotrigs ]
      else:
        notvetoed[col]   = list( trigs.getColumnByName( col ) )
        vetoed[col]      = list( vetotrigs.getColumnByName( col ) )
    except KeyError:
      err = 'Column %s not found in %s.' % ( col,triggers.tableName )
      raise KeyError, err

  # get z column
  if zcolumn:
    notvetoed[zcolumn] = list( trigs.getColumnByName( zcolumn ) )
    vetoed[zcolumn]    = list( vetotrigs.getColumnByName( zcolumn ) )

    if logz:
      notvetoed[zcolumn] = [math.log( t,10 ) for t in notvetoed[zcolumn]]
      vetoed[zcolumn]    = [math.log( t,10 ) for t in vetoed[zcolumn]]
  # or make a default replacement 
  else:
    notvetoed[zcolumn] = [1]*len( trigs )
    vetoed[zcolumn] = [1]*len( vetotrigs )

  # =============
  # generate plot
  # =============
  # fix flag for use with latex
  flag = flag.replace( '''_''','''\_''' )
  etg      = etg.replace( '''_''','''\_''' )

  # customise plot appearance
  pylab.rcParams.update( {"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 20,
                         "ytick.labelsize": 20,
                         "axes.titlesize": 16,
                         "axes.labelsize": 16,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20} )

  fig = pylab.figure( figsize=[12,6] )
  ax  = fig.gca()
  m = 5
  plots = []

  # define colour colormap
  cdict = matplotlib.cm.jet._segmentdata
  cmap = matplotlib.colors.LinearSegmentedColormap( 'clrs', cdict )
  # change start blue and end red to better colours
  cdict['red']=( [x == cdict['red'][-1] and ( 1,1,1 ) or x for x in cdict['red']] )
  cdict['blue']=( [x== cdict['blue'][0] and ( 0,1,1 ) or x for x in cdict['blue']] )

  # plot triggers not vetoed
  if segments is not None:  label = 'Not vetoed'
  else: label=None

  if len( notvetoed[xcolumn] )>=1:
    pylab.scatter( notvetoed[xcolumn], notvetoed[ycolumn],\
                   c=notvetoed[zcolumn], marker='o', cmap=cmap, label=label )

  #  construct colorbar
  if zcolumn:
    if len( notvetoed[zcolumn]+vetoed[zcolumn] ) >= 1:
      cmin = int( math.floor( min( notvetoed[zcolumn]+vetoed[zcolumn] ) ) )
      cmax = int( math.ceil( max( notvetoed[zcolumn]+vetoed[zcolumn] ) ) )
      colorticks = arange( math.floor( cmin ),math.ceil( cmax ),( cmax-cmin )/5 )
      # if colouring by SNR, move to standard DQ range of 5->100
      if zcolumn.lower()=='snr':
        cmin = math.log( 5,10 )
        cmax = math.log( 110,10 )
        colorticks = [1,1.5,2]
      # draw colorbar
      cb = pylab.colorbar( format = logz and "$10^{%.1f}$",ticks=colorticks )
      cb.set_clim( cmin,cmax )
      cb.draw_all()
      cb.ax.set_ylabel( zcolumn.replace( '_','\_' ).title() )

      # get loudest event and plot as gold star
      maxidx = list( notvetoed[zcolumn]+vetoed[zcolumn] )\
                    .index( max( notvetoed[zcolumn]+vetoed[zcolumn] ) )
      maxx   = list( notvetoed[xcolumn]+vetoed[xcolumn] )[maxidx]
      maxy   = list( notvetoed[ycolumn]+vetoed[ycolumn] )[maxidx]
      maxz   = math.pow( 10,list( notvetoed[zcolumn]+vetoed[zcolumn] )[maxidx] )

      pylab.plot( [maxx],[maxy],color='gold',marker='*',markersize=15 )

    else:
      maxx = 0
      maxy = 0
      maxz = 0

  # plot vetoed triggers if required
  if segments is not None:
    # if not vetoed triggers, plot something to ensure the label appears
    if len( vetoed[xcolumn] )<1:
      pylab.scatter( [notvetoed[xcolumn][0]],[notvetoed[ycolumn][0]],\
                    marker='x',label='Vetoed',visible=False )
    else:
      if zcolumn:
        c='k'
      else:
        c='r'
      pylab.scatter( vetoed[xcolumn],vetoed[ycolumn],marker='x',\
                     label='Vetoed',edgecolor=c )

  # set axes
  if logx:
    ax.set_xscale( 'log' )
  if logy:
    ax.set_yscale( 'log' )

  # print legend
  #if segments is not None:
  #  ax.legend()

  # get start time in UTC for axis
  if re.search( 'time',xcolumn+ycolumn ):
    zerostring = datetime( *date.XLALGPSToUTC( LIGOTimeGPS( zero ) )[:6] )\
                     .strftime( "%B %d %Y, %H:%M:%S %ZUTC" )
  # set x label and lim
  if re.search( 'time', xcolumn ):
    ax.set_xlabel( 'Time (%s) since %s (%s)'\
                  % ( t_string[t_unit], zerostring, zero ) )
    pstart = float( start-zero )/t_unit
    pend   = float( end-zero )/t_unit
    ax.set_xlim( pstart, pend)
  else:
    ax.set_xlabel( '%s' % ( xcolumn.replace( '_',' ' ).title() ) )
    if len( triggers )>=1:
      ax.set_xlim( min( vetoed[xcolumn]+notvetoed[xcolumn] )*0.99,\
                   max( vetoed[xcolumn]+notvetoed[xcolumn] )*1.01 )
  # set y label and lim
  if re.search( 'time',ycolumn ):
    ax.set_ylabel( 'Time (%s) since %s (%s)'\
                  % ( t_string[t_unit], zerostring, zero ) )
    if len( triggers )>=1:
      ax.set_ylim( 0,float( end-start )/t_unit )
  else:
    ax.set_ylabel( '%s' % ( ycolumn.replace( '_',' ' ).title() ) )
    if len( triggers )>=1:
      ax.set_ylim( min( vetoed[ycolumn]+notvetoed[ycolumn] )*0.99,\
                  max( vetoed[ycolumn]+notvetoed[ycolumn] )*1.01 )

  if ylim:
    ax.set_ylim( tuple(ylim) )

  # set title
  tit = '%s triggers' % ( etg )
  if segments is not None:
    tit += ' \&  %s segments' % ( flag )
  tit += ': %s-%s' % ( start, end )
  ax.set_title( tit, x=0.5, y=1.03 )

  # set subtitle
  if zcolumn:
    if re.search( 'time',xcolumn ):  maxx = maxx*t_unit+start
    if re.search( 'time',ycolumn ):  maxy = maxy*t_unit+start

    subtit = 'Loudest event: %s=%s %s=%.2f %s=%.2f'\
             % ( xcolumn.replace( '_','\_' ),maxx,\
                ycolumn.replace( '_','\_' ),maxy,\
                zcolumn.replace( '_','\_' ),maxz )
    ax.text( 0.5,1.001,subtit,horizontalalignment='center',\
            transform = ax.transAxes )

  # get both major and minor grid lines
  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )
  ax.set_axisbelow( True )
  fig.savefig( outfile )

# =============================================================================
# Plot science segment histogram
# =============================================================================

def plot_segment_hist( segments,outfile,flag=None,coltype=int,\
                      logx=False,logy=False ):

  """
    Plots a histogram of segment duration for the glue.segments.segmentlist
    segments.

    Arguments:

      segments : glue.segments.segmentlist
        list of segments with which to veto triggers
      outfile : string 
        string path for output plot
   
    Keyword arguments:

      flag : string
        display name for segments, normally the name of the DQ flag
      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
  """

  durations = [seg.__abs__() for seg in segments]

  # set numbins
  if not len( segments )<1:
    if coltype==int and not logx:
      numbins = min( max( durations )-min( durations ),200 )
    else:
      numbins = min( len( segments ),200 )

  if logx:
    durations = [math.log10( d ) for d in durations]

  # fix flag for use with latex
  if flag:
    flag = flag.replace( '''_''','''\_''' )

  # customise plot appearance
  pylab.rcParams.update( {"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16,
                         "axes.titlesize": 20,
                         "axes.labelsize": 16,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20} )

  # generate plot
  fig = pylab.figure( figsize=[12,6] )
  ax  = fig.gca()

  # calculate histogram
  if not len( segments )<1:

    hrange=( min( durations ),max( durations ) )

    n,b,p = pylab.hist( durations,bins=numbins,range=hrange,\
                       histtype='bar',visible=False )

    bins=[]
    barwidth=[]

    if logx:
      for i in range( len( b )-1 ):
        bins.append( math.pow( 10,b[i] ) )
        barwidth.append( math.pow( 10,b[i+1] )-math.pow( 10,b[i] ) )
    else:
      bins = b[0:-1]
      try:
        barwidth = [bins[1]-bins[0]]*len( bins )
      except IndexError:
        barwidth = [1] 
 
  else:
    bins = []
    n = []
    barwidth = []

  # fix base for logy
  if logy:
    base = 0.9
  else:
    base = 0
  
  n = [c-base for c in n]
  
  # plot histogram
  ax.bar( bins,n,width=barwidth,bottom=base,\
         color='blue',edgecolor='black',log=logy )

  # set axes
  if logx:
    ax.set_xscale( 'log' )
  if logy:
    ax.set_yscale( 'log' )

  if len( segments )>=1:
    if logx:
      # reset limits to logscale
      ax.set_xlim( ( math.pow( 10,min( durations ) )*0.99 ),\
                  ( math.pow( 10,max( durations ) )*1.01 ) )

    ax.set_ylim( base,math.pow( 10,math.log10( ( max( n )+base )*1.01 ) ) )
  ax.set_xlabel( 'Length of segment ( seconds )' )
  ax.set_ylabel( 'Number of segments' )
  tit = 'Segment Duration Histogram'
  if flag:
    tit = '%s %s' % ( flag,tit )
  ax.set_title( tit )

  # get both major and minor grid lines
  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )
  ax.set_axisbelow( True )
  fig.savefig( outfile )

