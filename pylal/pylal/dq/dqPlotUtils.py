#!/usr/bin/env python

from __future__ import division
import types,sys,os,math,numpy,subprocess,re

from glue.ligolw import lsctables
from glue.segments import segment, segmentlist
from glue import segmentsUtils

from datetime import datetime
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal import date

from matplotlib import mlab, use
use('Agg')
import pylab

from pylal.dq import dqDataUtils,dqSegmentUtils

# compile xml.gz search
xml = re.compile('(xml$|xml.gz$)')

# =============================================================================
# Function to execute shell command and get output
# =============================================================================
def make_external_call(cmd,shell=False):
  args = shlex.split(str(cmd))
  p = subprocess.Popen(args,shell=shell,\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p_out, p_err = p.communicate()
  if p.returncode != 0:
    raise ValueError, "Command %s failed. Stderr Output: \n%s" %( cmd, p_err)

  return p_out, p_err

# =============================================================================
# Function to plot before/after cumulative SNR histograms
# =============================================================================
def plot_veto_hist(start,end,trigs,segs,outfile,livetime=None,\
                   vetoname='unknown segments',\
                   column='snr',etg='unkown',fill=False,logx=True):

  # set up outfile
  outfile = os.path.abspath(outfile)
  outdir = os.path.split(outfile)[0]
  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  if not livetime:
    livetime=end-start

  # fix name for latex
  vetoname = vetoname.replace('_','\_')

  # customise plot appearance
  pylab.rcParams.update({"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 12,
                         "ytick.labelsize": 12,
                         "axes.titlesize": 20,
                         "axes.labelsize": 18,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20})

  # generate vetoed trigger list 
  aftertrigs  = trigs.veto(segs)

  # set up histogram data
  try:
    if logx:
      start_data = [math.log10(trig.__getattribute__(column))\
                    for trig in trigs]
      end_data   = [math.log10(trig.__getattribute__(column))\
                    for trig in aftertrigs]
    else:
      start_data = [trig.__getattribute__(column) for trig in trigs]
      end_data   = [trig.__getattribute__(column) for trig in aftertrigs]
  except KeyError:
    err = 'Column %s not found in %s.' % (column,trigs.tableName)
    raise KeyError, err

  # generate histogram
  num_bins=1000
  fig = pylab.figure(figsize=[12,6])
  ax  = fig.gca()
  ax.loglog()

  if len(trigs)>=1:
    start_n,start_bins,start_p = pylab.hist(start_data,bins=num_bins,\
                                            range=(min(start_data),\
                                            max(start_data)),\
                                            histtype='stepfilled',\
                                            cumulative=-1,facecolor='red',\
                                            visible=False)
    end_n,end_bins,end_p = pylab.hist(end_data,bins=num_bins,\
                                      range=(min(start_data),max(start_data)),\
                                      histtype='stepfilled',\
                                      cumulative=-1,facecolor='green',\
                                      visible=False)

  bins = []
  for i in range(len(start_bins)-1):
    if logx:
      bins.append(0.5*(math.pow(10,start_bins[i])+math.pow(10,start_bins[i+1])))
    else:
      bins.append(0.5*(start_bins[i]+start_bins[i+1]))

  # reset zero values to base
  livetime = float(livetime)
  base = 0.5/livetime
  def zero_replace(num,base):
    if num==0:
      return base
    else:
      return num
  start_n = map(lambda n: zero_replace(n,base), start_n)
  end_n   = map(lambda n: zero_replace(n,base), end_n)

  # convert number to rate
  start_n = [n/livetime for n in start_n]
  end_n   = [n/livetime for n in end_n]

  # plot data
  if not fill:
    ax.plot(bins,start_n,'r',linewidth=2,label='Before vetoes')
    ax.plot(bins,end_n,'g',linewidth=2,label='After vetoes')
  if fill:
    ax.plot(bins,start_n,'r',linewidth=0,label='Before vetoes')
    ax.plot(bins,end_n,'g',linewidth=0,label='After vetoes')
    ax.fill_between(bins,base,start_n,color='r',edgecolor='k',linewidth=0.5)
    ax.fill_between(bins,base,end_n,color='g',edgecolor='k',linewidth=0.5,\
                    alpha=0.9)

  # figure sundries
  vetoname = vetoname
  leg = ax.legend()
  for l in leg.get_lines():
    l.set_linewidth(4)
  ax.set_xlim(min(bins),max(bins))
  ax.set_ylim(base,max(start_n)*1.01)
  ax.set_xlabel('SNR')
  ax.set_ylabel('Cumulative rate (Hz)')
  ax.set_title('Effect of '+vetoname+' on '+etg+' triggers')
  ax.grid(True,which='major')
  ax.grid(True,which='majorminor')
  ax.set_axisbelow(True)
  fig.savefig(outfile)
  #ax.close()

# =============================================================================
# Function to plot SNR vs. time with veto segments
# =============================================================================
def plot_veto_trigs(start,end,trigs,segs,outfile,vetoname='unknown segments',\
                    xcolumn='time',ycolumn='snr',etg='unknown'):

  """
    Plots ycolumn against xcolumn for columns in given
    Sngl{Burst,Inspiral}Table object trigs, highlighting those entries falling
    inside one of the entries in the glue.segments.segmentlist object segs. 

    'time' given as a column name is a special case, since s and ns times are
    stored separately in the SnglTable structures. In this case the
    trigger.get_xxxx() function is called.

    Keyword arguments:
    
      start: float, numerical GPS time
      end: float, numerical GPS time
      cache: list, trigger file paths
      segfile: str, segment xml or segwizard file containing veto segments
      outfile: str, output filepath
      segdeffile: str, segment xml or segwizard file containing segment summary
      vetoname: str, display name for veto segments
      etg: ['ihope' | 'omega']
      xcolumn: str, any compatible column to plot on the x-axis
      ycolumn: str, any compatible column to plot on the y-axis
      table: str, name of table to extract from trigger xml files

  """

  # set up outfile
  outfile = os.path.abspath(outfile)
  outdir = os.path.split(outfile)[0]
  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  # fix vetoname for use with latex
  vetoname = vetoname.replace('''_''','''\_''')

  # customise plot appearance
  pylab.rcParams.update({"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 12,
                         "ytick.labelsize": 12,
                         "axes.titlesize": 20,
                         "axes.labelsize": 16,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20})

  # set up vetoed/nonvetoed trigger lists
  if re.search('sngl_inspiral',trigs.tableName):
    trigs_not  = trigs.veto(segs)
    trigs_veto = trigs.vetoed(segs)
  else:
    trigs_not = [t for t in trigs if t.get_peak() not in segs]
    trig_veto = [t for t in trigs if t.get_peak() in segs]

  # set plot time unit
  if (end-start) < 20000:
    t_unit = 60
  elif (end-start) >= 20000 and (end-start) < 604800:
    t_unit = 3600
  else:
    t_unit = 86400
  t_string = {60:'minutes',3600:'hours',86400:'days'}

  # set up plot lists
  try:
    if xcolumn=='time':
      if re.search('sngl_inspiral',trigs.tableName):
        x_not  = [float(trig.get_end()-start)/t_unit for trig in trigs_not]
        x_veto = [float(trig.get_end()-start)/t_unit for trig in trigs_veto]
      else:
        x_not  = [float(trig.get_peak()-start)/t_unit for trig in trigs_not]
        x_veto = [float(trig.get_peak()-start)/t_unit for trig in trigs_veto]
    else:
      x_not  = trigs_not.get_column(xcolumn).tolist()
      x_veto = trigs_veto.get_column(xcolumn).tolist()
  except KeyError:
    err = 'Column %s not found in %s.' % (xcolumn,trigs.tableName)
    raise KeyError, err

  try:
    if ycolumn=='time':
      if re.search('sngl_inspiral',trigs.tableName):
        y_not  = [float(trig.get_end()-start)/t_unit for trig in trigs_not]
        y_veto = [float(trig.get_end()-start)/t_unit for trig in trigs_veto]
      else:
        y_not  = [float(trig.get_peak()-start)/t_unit for trig in trigs_not]
        y_veto = [float(trig.get_peak()-start)/t_unit for trig in trigs_veto]
    else:
      y_not  = trigs_not.get_column(ycolumn).tolist()
      y_veto = trigs_veto.get_column(ycolumn).tolist()
  except KeyError:
    err = 'Column \'%s\' not found in \'%s.\'' % (ycolumn,trigs.tableName)
    raise KeyError, err

  
  # generate plot
  fig = pylab.figure(figsize=[12,6])
  ax  = fig.gca()
  # plot x versus y in legend markersize
  m = 5
  p1 = ax.plot(x_not,y_not,'b.',label='Not vetoed',\
               markersize=3*m)
  p2 = ax.plot(x_veto,y_veto,'r.',label='Vetoed',\
               markersize=3*m)

  # rescale plot markers down to normal size
  for p in [p1,p2]:
    p = p[0]
    p.set_markersize(m)

  # sundries
  ax.semilogy()
  ax.legend()
  if re.search('time',xcolumn+ycolumn):
    startstring = datetime(*date.XLALGPSToUTC(start)[:6])\
                .strftime("%B %d %Y, %H:%M:%S %ZUTC")
  if re.search('time',xcolumn):
    ax.set_xlabel('Time ('+t_string[t_unit]+') since '+str(startstring))
    ax.set_xlim(0,float(end-start)/t_unit)
  else:
    ax.set_xlabel('%s' % (xcolumn.replace('_',' ').title()))
    if len(trigs)>=1:
      ax.set_xlim(min(x_veto+x_not)*0.99,max(x_veto+x_not)*1.01)
  if re.search('time',ycolumn):
    ax.set_ylabel('Time ('+t_string[t_unit]+') since '+str(startstring))
    if len(trigs)>=1:
      ax.set_ylim(0,float(end-start)/t_unit)
  else:
    ax.set_ylabel('%s' % (ycolumn.replace('_',' ').title()))
    if len(trigs)>=1:
      ax.set_ylim(min(y_veto+y_not)*0.99,max(y_veto+y_not)*1.01)
  ax.set_title('Effect of '+vetoname+' on '+etg+' triggers')

  # get both major and minor grid lines
  ax.grid(True,which='major')
  ax.grid(True,which='majorminor')
  ax.set_axisbelow(True)
  fig.savefig(outfile)

# =============================================================================
# Plot science segment histogram
# =============================================================================
def plot_segment_hist(segs,outfile,name=None,logx=False,logy=False):

  """Plots a histogram of segment duration for the given glue.segmentlist object segs."""

  # set up outfile
  outfile = os.path.abspath(outfile)
  outdir = os.path.split(outfile)[0]
  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  durations = [seg.__abs__() for seg in segs]
  # set times
  t_unit = 1
  t_string = {1:'seconds',60:'minutes',3600:'hours',86400:'days'}

  durations = [float(d)/t_unit for d in durations]
  if logx:
    durations = [math.log10(d) for d in durations]

  # fix vetoname for use with latex
  if name:
    name = name.replace('''_''','''\_''')

  # customise plot appearance
  pylab.rcParams.update({"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 5,
                         "xtick.labelsize": 12,
                         "ytick.labelsize": 12,
                         "axes.titlesize": 20,
                         "axes.labelsize": 16,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 20})

  # generate plot
  fig = pylab.figure(figsize=[12,6])
  ax  = fig.gca()

  # calculate histogram
  if not len(segs)<1: 
    numbins = min(len(segs),200)
    n,b,p = pylab.hist(durations,bins=numbins,\
                       range=(min(durations),max(durations)),\
                       histtype='bar',visible=False)
  
    bins=[]
    barwidth=[]
    if logx:
      for i in range(len(b)-1):
        bins.append(math.pow(10,b[i]))
        barwidth.append(math.pow(10,b[i+1])-math.pow(10,b[i]))
    else:
      bins = b[0:-1]
      barwidth = [bins[1]-bins[0]]*len(bins)
      
 
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
  ax.bar(bins,n,width=barwidth,bottom=base,\
         color='blue',edgecolor='black',log=logy)

  if len(segs)>=1:
    if logx:
      ax.semilogx()
      ax.set_xlim(min(bins)*0.99,(max(bins)+max(barwidth))*1.01)
    else:
      ax.set_xlim(0,(max(bins)+max(barwidth))*1.01)
    ax.set_ylim(base,math.pow(10,math.log10((max(n)+base)*1.01)))
  ax.set_xlabel('Length of segment (%s)' %(t_string[t_unit]))
  ax.set_ylabel('Number of segments')
  tit = 'Segment Duration Histogram'
  if name:
    tit = '%s %s' % (name,tit)
  ax.set_title(tit)

  # get both major and minor grid lines
  ax.grid(True,which='major')
  ax.grid(True,which='majorminor')
  ax.set_axisbelow(True)
  fig.savefig(outfile)

