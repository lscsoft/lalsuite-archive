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
# Set plot parameters
# =============================================================================

def set_rcParams():

  # customise plot appearance
  pylab.rcParams.update( {"text.usetex": True,
                          "text.verticalalignment": "center",
                          "lines.linewidth": 5,
                          "xtick.labelsize": 18,
                          "ytick.labelsize": 18,
                          "axes.titlesize": 24,
                          "axes.labelsize": 20,
                          "axes.linewidth": 1,
                          "grid.linewidth": 1,
                          "legend.fontsize": 18,
                          "legend.loc": "best",
                          "figure.figsize": [12,6],
                          "figure.dpi": 80,
                          "axes.grid": True,
                          "axes.axisbelow": True } )

# =============================================================================
# Plot before/after cumulative SNR histograms
# =============================================================================

def plot_trigger_hist( triggers, outfile, column='snr', segments=None,\
                       start=None, end=None, bins=1000,\
                       flag='unknowns', etg='Unknown',\
                       livetime=None, fill=False, logx=True, logy=True,
                       cumulative=True, rate=True, xlim=None, greyscale=False ):

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

    if not start and not end:
      times = [ get_time( t ) for t in triggers ] 
    if not start:
      start = min( times )
    if not end:
      end   = max( times )
    livetime = end-start 
  livetime = float(livetime)

  # generate vetoed trigger list: inspiral trigs have dedicated function 
  if not segments:
    segments = segmentlist()

  aftertriggers = table.new_from_template( triggers )
  aftertriggers.extend([ t for t in triggers if get_time(t) not in segments ])

  afterlivetime = livetime - float(segments.__abs__())

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

  if segments:
    color = ['r','g']
  else:
    color = ['b']
  if greyscale:
    color = ['k','k']
    linestyle = ['-','--']
  else:
    linestyle = ['-','-']

  for i,c in enumerate(color):
    color[i] = '%s%s' % ( c, linestyle[i] )

  # generate histogram
  if not bins:
    bins = 1000 

  if len( triggers )>=1:

    if cumulative:
      cumulative = -1

    start_n,start_bins,start_p = pylab.hist( start_data, bins=bins,\
                                             range=( min( start_data ),\
                                                     max( start_data ) ),\
                                             histtype='stepfilled',\
                                             cumulative=cumulative,\
                                             visible=False )

  else:
    bins = []
    start_n = []
    
  if len( aftertriggers)>=1:
    end_n,end_bins,end_p = pylab.hist( end_data, bins=bins,\
                                       range=( min( start_data ),\
                                               max( start_data ) ),\
                                       histtype='stepfilled',\
                                       cumulative=cumulative,\
                                       visible=False )

  else:
    end_n = []

  if len( triggers )>=1:
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
      end_n   = [n/afterlivetime for n in end_n]

  # fix names for latex
  if flag:
    flag = flag.replace( '_','\_' )
  else:
    flag = 'unknown'
  column = ' '.join([ w.title().replace('Snr','SNR')\
                      for w in column.split('_') ])

  # customise plot appearance
  set_rcParams()

  # plot data
  fig = pylab.figure( figsize=[12,6] )
  ax  = fig.gca()
  if logx:
    ax.set_xscale('log')
  if logy:
    ax.set_yscale('log')

  if not fill:
    ax.plot( bins,start_n,color[0],linewidth=2,label='Before vetoes' )
    if len( aftertriggers )>=1:
      ax.plot( bins,end_n,color[-1],linewidth=2,label='After vetoes' )
  if fill:
    ax.plot( bins,start_n,color[0],linewidth=0,label='Before vetoes' )
    ax.fill_between( bins,base,start_n,color=color[-1][0],edgecolor='k',\
                     linewidth=0.5 )
    if len( aftertriggers )>=1:
      ax.plot( bins,end_n,color[-1],linewidth=0,label='After vetoes' )
      ax.fill_between( bins,base,end_n,color=color[-1][0],edgecolor='k',\
                       linewidth=0.5, alpha=0.9 )

  # figure sundries
  if segments:
    leg = ax.legend( loc='best' )
    for l in leg.get_lines():
      l.set_linewidth( 4 )
  if bins:
    ax.set_xlim( float(min(bins)), float(max(bins)) )
    ax.set_ylim( base,max( start_n )*1.01 )
  if xlim:
    ax.set_xlim( tuple(xlim) )
  ax.set_xlabel( column.replace( '_','\_' ) )
  if rate and cumulative:
    ax.set_ylabel( 'Cumulative rate (Hz)' )
  elif rate:
    ax.set_ylabel( 'Rate (Hz)' )
  elif not rate and cumulative:
    ax.set_ylabel( 'Cumulative number' )
  elif not rate and not cumulative:
    ax.set_ylabel( 'Number' )

  tit = '%s triggers' % ( etg.replace( '_','\_' ) )
  if segments:
    tit += ' and %s segments' % ( flag )
  ax.set_title( tit, x=0.5, y=1.035 )

  if start and end:
    subtit = '%s-%s' % ( start, end )
    ax.text( 0.5, 1.03, subtit, horizontalalignment='center',\
             transform = ax.transAxes, verticalalignment='top' )

  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )
  fig.savefig( outfile, bbox_inches='tight' )

# =============================================================================
# Function to plot SNR vs. time with veto segments
# =============================================================================

def plot_triggers( triggers, outfile, etg='Unknown',\
                   start=None, end=None, zero=None,\
                   segments=None, flag=None,\
                   xcolumn='time', ycolumn='snr', zcolumn=None,\
                   xlabel=None, ylabel=None, zlabel=None,\
                   logx=False, logy=True, logz=True, xlim=None, ylim=None,\
                   greyscale=False, set_plot_params=True ):

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

  span = segment([start,end])

  columns = [ xcolumn, ycolumn]
  xcolumn = xcolumn.lower()
  ycolumn = ycolumn.lower()
  if zcolumn:
    columns.append(zcolumn)
    zcolumn = zcolumn.lower()

  for i,c in enumerate(columns):
    columns[i] = ' '.join([ w.title().replace('Snr','SNR')\
                            for w in c.split('_') ])

  # sort triggers by z param
  if zcolumn:
    triggers.sort( key=lambda trig: trig.__getattribute__( zcolumn ),\
                   reverse=False )

  # apply veto segments if required
  if segments is not None:
    segments = segmentlist( segments )
    segments & [start,end]
  
    # set up vetoed/nonvetoed trigger lists, inspiral triggers have dicated func
    trigs =  table.new_from_template( triggers )
    vetotrigs = table.new_from_template( triggers )

    trigs.extend( t for t in triggers if get_time( t ) not in segments )
    vetotrigs.extend( t for t in triggers if get_time( t ) in segments )

  # or copy lists with empty vetoed triggers list
  else:
    trigs     = table.new_from_template( triggers)
    trigs.extend( t for t in triggers if get_time(t) in span )
    vetotrigs = table.new_from_template( triggers )

  # set plot time unit whether it's used or not
  if (end-start) < 1000:
    unit = 1
  elif (end-start) < 20000:
    unit = 60
  elif (end-start) >= 20000 and (end-start) < 604800:
    unit = 3600
  else:
    unit = 86400
  unitstr = {1:'seconds', 60:'minutes',3600:'hours',86400:'days'}

  if ( end-start ) < 20000:
    unit = 60
  elif ( end-start ) >= 20000 and ( end-start ) < 604800:
    unit = 3600
  else:
    unit = 86400
  unitstr = {60:'minutes',3600:'hours',86400:'days'}

  # set up plot lists
  notvetoed = {}
  vetoed = {}

  for col in [xcolumn,ycolumn]:
    try:
      # treat 'time as special case'
      if col=='time':
        notvetoed[col] = [ float( get_time(t) - zero ) / unit\
                           for t in trigs ]
        vetoed[col] = [ float( get_time(t) - zero ) / unit\
                        for t in vetotrigs ]
      else:
        notvetoed[col]   = list( trigs.getColumnByName( col ) )
        vetoed[col]      = list( vetotrigs.getColumnByName( col ) )
    except KeyError:
      err = 'Column %s not found in %s.' % ( col,triggers.tableName )
      raise KeyError, err

  # get z column
  if zcolumn:
    if logz:
      notvetoed[zcolumn] = [ math.log10(t) for t in\
                             trigs.getColumnByName( zcolumn ) ]
      vetoed[zcolumn]    = [ math.log10(t) for t in\
                             vetotrigs.getColumnByName( zcolumn ) ]
    else:
      notvetoed[zcolumn] = list( trigs.getColumnByName( zcolumn ) )
      vetoed[zcolumn]    = list( vetotrigs.getColumnByName( zcolumn ) )

  # or make a default replacement 
  else:
    notvetoed[zcolumn] = [1]*len( trigs )
    vetoed[zcolumn] = [1]*len( vetotrigs )

  numtrigs = len( trigs )
  numveto  = len( vetotrigs )

  # =============
  # generate plot
  # =============
  # fix flag for use with latex
  if flag:
    flag = flag.replace( '_', '\_' )
  else:
    flag = 'unknown'
  if etg:
    etg  = etg.replace( '_', '\_' )
  else:
    etg  = 'Unknown'

  # customise plot appearance
  if set_plot_params:
    set_rcParams()

  fig = pylab.figure( figsize=[12,6] )
  ax  = fig.gca()
  #ax  = pylab.axes([0.1, 0.1, 0.9, 0.8])
  m = 5
  plots = []

  # define colour colormap
  if greyscale:
    cdict = matplotlib.cm.hot._segmentdata
  else:
    cdict2 = matplotlib.cm.jet._segmentdata
    cdict = cdict2
    #cdict = {}
    #for c in cdict2.keys():
    #  row = map(list, cdict2[c][1:])
    #  row[0][0] = 0
    #  cdict[c] = map(tuple,row)
  cmap = matplotlib.colors.LinearSegmentedColormap( 'clrs', cdict )

  # define colour range
  cmin = None
  if zcolumn:
    if numtrigs + numveto >= 1:
      cmin = min( notvetoed[zcolumn]+vetoed[zcolumn] )
      cmax = max( notvetoed[zcolumn]+vetoed[zcolumn] )
      colorticks = arange( cmin, cmax, float(cmax - cmin)/5 )
    # if colouring by SNR, move to standard DQ range of 5->100
    if zcolumn.lower()=='snr' or numtrigs + numveto < 1:
      if logz:
        cmin = math.log( 3, 10 )
        cmax = math.log( 110, 10 )
        colorticks = [ math.log10(3), math.log10(10), math.log10(30),\
                       math.log10(100) ]
  if cmin==None:
    cmin = 3
    cmax = 100
    colorticks = [10,50,100]

  if numtrigs >= 1:

    p1 = ax.scatter( notvetoed[xcolumn], notvetoed[ycolumn],\
                     c=notvetoed[zcolumn], marker='o', cmap=cmap,\
                     vmin=cmin, vmax=cmax )

  # plot vetoed triggers if required
  if numveto >= 1:
    if zcolumn:
      c='k'
    else:
      c='r'
    p2 = ax.scatter( vetoed[xcolumn],vetoed[ycolumn], marker='x',\
                     label='Vetoed',edgecolor=c, vmin=cmin, vmax=cmax )

  if numtrigs < 1:
    p1 = ax.scatter( [0], [0], c=[cmin], cmap=cmap, vmin=cmin, vmax=cmax,\
                     visible=False )

  #  construct colorbar
  if zcolumn:

    if numtrigs + numveto >= 1:
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


    # draw colorbar
    formatter = matplotlib.ticker.FuncFormatter( lambda x,pos:\
                                                 "%d" % round(math.pow(10,x)) )
    cb = ax.figure.colorbar( p1, format = logz and formatter,\
                             ticks=colorticks )
    cb.draw_all()
    if not zlabel:
      zlabel = columns[2]
    cb.ax.set_ylabel( zlabel )

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
    if not xlabel:
      xlabel = 'Time (%s) since %s (%s)'\
                % ( unitstr[unit], zerostring, zero )
    pstart = float( start-zero )/unit
    pend   = float( end-zero )/unit
    ax.set_xlim( pstart, pend)

  else:
    if not xlabel:
      xlabel = columns[0]
    if len( triggers )>=1:
      ax.set_xlim( min( vetoed[xcolumn]+notvetoed[xcolumn] )*0.99,\
                   max( vetoed[xcolumn]+notvetoed[xcolumn] )*1.01 )

  ax.set_xlabel( xlabel )

  # set y label and lim
  if re.search( 'time',ycolumn ):
    if not ylabel:
      ylabel = 'Time (%s) since %s (%s)'\
               % ( unitstr[unit], zerostring, zero )
    if len( triggers )>=1:
      ax.set_ylim( 0,float( end-start )/unit )
  else:
    if not ylabel:
      ylabel = columns[1]
    if len( triggers )>=1:
      ax.set_ylim( min( vetoed[ycolumn]+notvetoed[ycolumn] )*0.99,\
                  max( vetoed[ycolumn]+notvetoed[ycolumn] )*1.01 )

  ax.set_ylabel( ylabel )

  if xlim:
    ax.set_xlim( tuple(xlim) )
  if ylim:
    ax.set_ylim( tuple(ylim) )

  # set title
  tit = '%s triggers' % ( etg )
  if segments is not None:
    tit += ' \&  %s segments' % ( flag )
  ax.set_title( tit, x=0.5, y=1.035 )

  # set subtitle
  if zcolumn and zcolumn.lower() == 'snr':
    if re.search( 'time',xcolumn ):  maxx = maxx*unit+start
    if re.search( 'time',ycolumn ):  maxy = maxy*unit+start

    subtit = 'Loudest event: %s=%s %s=%.2f %s=%.2f'\
             % ( columns[0], maxx,\
                 columns[1], maxy,\
                 columns[2], maxz )
  else:
    subtit = '%s-%s' % ( start, end )
  ax.text( 0.5, 1.03, subtit, horizontalalignment='center',\
           transform = ax.transAxes, verticalalignment='top' )

  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )

  # get both major and minor grid lines
  fig.savefig( outfile, bbox_inches='tight' )

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

  # fix flag and columns for use with latex
  if flag:
    flag = flag.replace( '_','\_' )

  # customise plot appearance
  set_rcParams()

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

  fig.savefig( outfile, bbox_inches='tight' )

# =============================================================================
# Rate versus time
# =============================================================================

def plot_trigger_rate( triggers, outfile, average=600, start=None, end=None,\
                       bincolumn='peak_frequency', bins=[], etg='Unknown',\
                       logx=False, logy=True, zero=None, ylim=None ):

  """
    INSERT DOCSTRING
  """

  get_time = def_get_time( triggers.tableName )

  # set start and end times
  if not start and not end:
    times = [ get_time( t ) for t in triggers ]
  if not start:
    start = min( times )
  if not end:
    end   = max( times )

  if not zero:
    zero = start

  # format bins
  bins = [ map( float, bin) for bin in bins ]

  bintrigs = {}
  label = {}
  if not bins:
    bins = [ (0,float('inf')) ]
    label[bins[0]] = None
  else:
    for i,bin in enumerate(bins):
      label[i] = '-'.join( map( str, bin ) )

  # get triggers for each bin 
  for i,bin in enumerate(bins):
    bintrigs[i] = table.new_from_template(triggers)
    bintrigs[i].extend([ t for t in triggers if\
                           bin[0] <= t.__getattribute__(bincolumn) < bin[1] ])

  # calculate rates
  xcol = []
  rate = {}
  for i,bin in enumerate(bins):
    rate[i] = []
  s = start
  while s < end:
    e = min( s+average, end )
    xcol.append(float(e+s)/2)
    for i,bin in enumerate(bins):
      numtrigs = len([ t for t in bintrigs[i] if s<=get_time(t)<e ])
      rate[i].append(numtrigs/average)
    s = e

  # customise plot appearance
  set_rcParams()

  # set plot time unit whether it's used or not
  if ( end-start ) < 20000:
    unit = 60
  elif ( end-start ) >= 20000 and ( end-start ) < 604800:
    unit = 3600
  else:
    unit = 86400
  unitstr = {60:'minutes',3600:'hours',86400:'days'}

  # set time axis
  xcol = [ float(t-zero)/unit for t in xcol ]

  # generate figure
  fig = pylab.figure( figsize=[12,6] )
  ax = fig.gca()

  # plot rates
  for i in range(len(bins)):
    ax.plot( xcol, rate[i], 'o', label=label[i] )

  if len(triggers)<1:
    ax.plot([1],[0.1],visible=False)

  # set axes
  if logx:
    ax.set_xscale( 'log' )
  if logy:
    ax.set_yscale( 'log' )

  # get start time in UTC for axis
  zerostring = datetime( *date.XLALGPSToUTC( LIGOTimeGPS( zero ) )[:6] )\
                   .strftime( "%B %d %Y, %H:%M:%S %ZUTC" )

  # plot sundries
  ax.legend(loc='best')
  ax.set_ylabel( 'Rate (Hz)' )
  ax.set_xlabel( 'Time (%s) since %s (%s)'#
                 % ( unitstr[unit], zerostring, zero ) )

  ax.set_xlim( float(start-zero)/unit, float(end-zero)/unit )
  if ylim:
    ax.set_ylim( tuple(ylim) )
  etg = etg.replace('_','\_')
  bincolumn = ' '.join([ w.title().replace('Snr','SNR')\
                         for w in bincolumn.split('_') ])


  tit = '%s triggers' % (etg)
  if label[0]:
    tit += ' binned by %s' % bincolumn
  ax.set_title( tit, x=0.5, y=1.035 )
  subtit = '%s-%s' % ( start, end )
  ax.text( 0.5, 1.03, subtit, horizontalalignment='center',
           transform=ax.transAxes, verticalalignment='top' )

  fig.savefig( outfile, bbox_inches='tight' )

# =============================================================================
# Plot time series
# =============================================================================

def plot_time_series( data, outfile, start=None, end=None, zero=None, \
                      zeroindicator=False, subtitle=None, style='-',\
                      logx=False, logy=False, ylim=None, ylabel=None ):

  """
    Plot the time series of a given set (or given sets) of data.

    Arguments:

      data : list
        list of (ChannelName,time,amplitude) tuples with channel name (or data 
        source) and time/amplitude arrays for each channel. Channels are
        plotted in the order given.
      outfile : str
        output plot path

    Keyword Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of plot
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of plot
      zero : [ float | int | LIGOTimeGPS ]
        time around which to centre plot
      zeroindicator : [ False | True ]
        indicate zero time with veritcal dashed line, default: False
      ylabel : str
        title string for y axis, default: 'Signal'
      subtitle : str
        descriptive string used for sub title
      style : str
        matplotlib style string for line formatting, '-','--','.','x', etc,
        default: '-'
  """

  # format times
  if not start:
    start = min([ min(d[1]) for d in data ])
  if not end:
    end   = max([ max(d[1]) for d in data ])
  if not zero:
    zero = start

  start = LIGOTimeGPS( start )
  end   = LIGOTimeGPS( end )
  zero  = LIGOTimeGPS( zero )


  # set plot time unit whether it's used or not
  if (end-start) < 1000:
    unit = 1
  elif (end-start) < 20000:
    unit = 60
  elif (end-start) >= 20000 and (end-start) < 604800:
    unit = 3600
  else:
    unit = 86400
  unitstr = {1:'seconds', 60:'minutes',3600:'hours',86400:'days'}

  # customise plot appearance
  set_rcParams()

  # plot data
  fig = pylab.figure()
  ax  = fig.gca()

  for channel,time,amplitude in data:

    time = [ float( t - zero ) / unit for t in time ]
    lab = str(channel).replace('_','\_')

    if style in ['-','--']:
      ax.plot( time, amplitude, style, label = lab,\
               linewidth = 0.7 )

    else:
      ax.plot( time, amplitude, style, label = lab,\
               markersize = 2 )

  ax.legend( loc = 'best' )

  if ylim:
    ax.set_ylim( tuple(ylim) )

  # FIXME add zero indicator
  axis_lims = ax.get_ylim()
  if zeroindicator:
    ax.plot([0,0],[axis_lims[0],axis_lims[1]],'r--', linewidth=2.0)
    ax.set_ylim([ axis_lims[0], axis_lims[1] ])

  # set x axis
  ax.set_xlim([ float(start-zero)/unit, float(end-zero)/unit ])
  zero = LIGOTimeGPS('%.3f' % zero )
  if zero.nanoseconds==0:
    zerostr = datetime( *date.XLALGPSToUTC(LIGOTimeGPS(zero ))[:6])\
                  .strftime( "%B %d %Y, %H:%M:%S %ZUTC" )
  else:
    zerostr = datetime( *date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                  .strftime( "%B %d %Y, %H:%M:%S %ZUTC" )
    zerostr = zerostr.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
  ax.set_xlabel( 'Time (%s) since %s (%s)' % ( unitstr[unit], zerostr, zero ) )

  # set y axis
  if ylabel:
    ax.set_ylabel(ylabel)
  else:
    ax.set_ylabel('Signal')

  # set legend
  leg = ax.legend()
  for l in leg.get_lines():
    l.set_linewidth( 4 )

  # set scale
  if logx:
    ax.set_xscale('log')
  if logy:
    ax.set_yscale('log')

  # set title
  ax.set_title( 'Time series', x=0.5, y=1.035 )
  if subtitle:
    ax.text( 0.5, 1.03, subtitle, horizontalalignment='center',
             transform=ax.transAxes, verticalalignment='top' )

  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )

  fig.savefig( outfile, bbox_inches = 'tight' )

# =============================================================================
# Plot spectrum
# =============================================================================

def plot_spectrum( data, outfile, logx=False, logy=False, xlim=None, ylim=None,\
                   subtitle=None, style='-', ylabel=None ):

  """
    Plot the time series of a given set (or given sets) of data.

    Arguments:

      data : list
        list of (channel,frequency,spectral amplitude) tuples with channel 
        name (or data source) and frequency/spectral amplitude arrays for each
        channel. Channels are plotted in given order.
      outfile : str
        output plot path

    Keyword Arguments:

      logx : [ False | True ]
        plot x axis in log scale, default False
      logy : [ False | True ]
        plot y axis in log scale, default False
      xlim : [ list | tuple ]
        (min,max) tuple/list with limits for x axis
      ylim : [ list | tuple ]
        (min,maxy) tuple/list with limits for y axis
      subtitle : str
        descriptive string used for sub title
      style : str
        matplotlib style string for line formatting, '-','--','.','x', etc,
        default: '-'
  """

  # customise plot appearance
  set_rcParams()

  # plot data
  fig = pylab.figure()
  ax  = fig.gca()

  for channel,freq,spec in data:

    lab = str(channel).replace('_','\_')

    if style in ['-','--']:
      ax.plot( freq, spec, style, label = lab,\
               linewidth = 0.7 )

    else:
      ax.plot( freq, spec, style, label = lab,\
               markersize = 2 )

  # set legend
  leg = ax.legend()
  for l in leg.get_lines():
    l.set_linewidth( 4 )

  # set scale
  if logx:
    ax.set_xscale('log')
  if logy:
    ax.set_yscale('log')

  # set limits
  if xlim:
    ax.set_xlim( tuple(xlim) )
  if ylim:
    ax.set_ylim( tuple(ylim) )

  # set labels
  ax.set_xlabel( 'Frequency (Hz)' )
  if ylabel:
    ax.set_ylabel(ylabel)
  else:
    ax.set_ylabel( 'Spectrum' )

  # set title
  ax.set_title( 'Frequency Spectrum', x=0.5, y=1.035 )
  if subtitle:
    ax.text( 0.5, 1.03, subtitle, horizontalalignment='center',
             transform=ax.transAxes, verticalalignment='top' )

  ax.grid( True,which='major' )
  ax.grid( True,which='majorminor' )

  fig.savefig( outfile, bbox_inches = 'tight' )

