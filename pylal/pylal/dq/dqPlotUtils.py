#!/usr/bin/env python

from __future__ import division
import types,sys,os,math,numpy,subprocess
import os

from glue.segments import segment, segmentlist
from glue import segmentsUtils

from datetime import datetime
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal import date

from matplotlib import mlab, use
use('Agg')
import pylab
from pylal.dq import dqDataUtils,dqSegmentUtils

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
def plot_veto_hist(start,end,cache,segfile,outfile,\
                   segdeffile=None,livetime=None,vetoname='VETO',\
                   etg='ihope',logx=True,table='sngl_inspiral',fill=False,\
                   verbose=False):

  #== set up outfile
  outfile = os.path.abspath(outfile)
  outdir = os.path.split(outfile)[0]
  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  #== set up params
  if isinstance(segfile,str):
    segfile = segfile.split(',')
  if isinstance(segdeffile,str):
    segdeffile = segdeffile.split(',')
  if not livetime:
    livetime=end-start

  #== customise plot appearance
  pylab.rcParams.update({
                         #"font.family": "serif",
                         "text.usetex": True,
                         "text.verticalalignment": "center",
                         #"lines.markersize": 12,
                         #"lines.markeredgewidth": 2,
                         "lines.linewidth": 5,
                         "font.size": 20,
                         "axes.titlesize": 20,
                         "axes.labelsize": 18,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16,
                         "legend.fontsize": 20})

  #== load triggers
  if verbose:
    print >>sys.stdout, "Loading triggers..."
  if isinstance(cache,str):
    cache = [cache]
  beforetrigs = []
  for f in cache:
    if f.endswith('xml'):
      beforetrigs += dqDataUtils.fromtrigxml(open(f,'r'),start,end,table)
    else:
      beforetrigs += dqDataUtils.fromtrigfile(open(f,'r'),start,end,etg)

  #== load segments
  if verbose:
    print >>sys.stdout, "Loading segments..."
  segs = segmentlist()
  if segfile:
    for f in segfile:
      f = os.path.realpath(f)
      #== check file exists
      if not os.path.isfile(f):
        print >>sys.stderr, "Error: Cannot find segment file "+f
        continue

      #== load segment definer file
      if segdeffile:
        idx = segfile.index(f)
        try:
          if segdeffile[idx].endswith('xml'):
            segdefs = dqSegmentUtils.fromsegmentxml(open(segdeffile[idx],'r'))
          else:
            segdefs = segmentsUtils.fromsegwizard(open(segdeffile[idx],'r'))
        except:
          segdefs = segmentlist([segment(start,end)])
      else:
        segdefs = segmentlist([segment(start,end)])
      if f.endswith('xml'):
        segs += dqSegmentUtils.fromsegmentxml(open(f,'r'))&segdefs
      else:
        segs += segmentsUtils.fromsegwizard(open(f,'r'))&segdefs
  if verbose:
    print >>sys.stdout, "Parsing vetoed triggers..."
  aftertrigs = dqDataUtils.parse_trigs(beforetrigs,segs,inclusive=False)
  #== set up histogram data
  if logx:
    start_snrs = [math.log10(trig.snr) for trig in beforetrigs]
    end_snrs   = [math.log10(trig.snr) for trig in aftertrigs]
  else:
    start_snrs = [trig.snr for trig in beforetrigs]
    end_snrs   = [trig.snr for trig in aftertrigs]

  #== generate histogram
  if verbose:
    print >>sys.stdout, "Generating histogram..."
  num_bins=1000
  fig = pylab.figure(figsize=[12,6])
  ax  = fig.gca()
  ax.loglog()

  start_n,start_bins,start_p = pylab.hist(start_snrs,bins=num_bins,\
                                          range=(min(start_snrs),\
                                          max(start_snrs)),\
                                          histtype='stepfilled',\
                                          cumulative=-1,facecolor='red',\
                                          visible=False)
  end_n,end_bins,end_p = pylab.hist(end_snrs,bins=num_bins,\
                                    range=(min(start_snrs),max(start_snrs)),\
                                    histtype='stepfilled',\
                                    cumulative=-1,facecolor='green',\
                                    visible=False)

  bins = []
  for i in range(len(start_bins)-1):
    if logx:
      bins.append(0.5*(math.pow(10,start_bins[i])+math.pow(10,start_bins[i+1])))
    else:
      bins.append(0.5*(start_bins[i]+start_bins[i+1]))

  #== reset zero values to base
  livetime = float(livetime)
  base = 0.5/livetime
  def zero_replace(num,base):
    if num==0:
      return base
    else:
      return num
  start_n = map(lambda n: zero_replace(n,base), start_n)
  end_n   = map(lambda n: zero_replace(n,base), end_n)

  #== convert number to rate
  start_n = [n/livetime for n in start_n]
  end_n   = [n/livetime for n in end_n]
  if verbose:
    print >>sys.stdout, "Plotting..."
  if not fill:
    ax.plot(bins,start_n,'r',linewidth=2,label='Before vetoes')
    ax.plot(bins,end_n,'g',linewidth=2,label='After vetoes')
  if fill:
    ax.plot(bins,start_n,'r',linewidth=0,label='Before vetoes')
    ax.plot(bins,end_n,'g',linewidth=0,label='After vetoes')
    ax.fill_between(bins,base,start_n,color='r',edgecolor='k',linewidth=0.5)
    ax.fill_between(bins,base,end_n,color='g',edgecolor='k',linewidth=0.5)

  #== figure sundries
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
  fig.savefig(outfile)
  #ax.close()
  if verbose:
    print >>sys.stdout, "Done."

# =============================================================================
# Function to plot SNR vs. time with veto segments
# =============================================================================
def plot_veto_trigs(start,end,cache,segfile,outfile,\
                    segdeffile=None,vetoname='VETO',etg='ihope',\
                    table='sngl_inspiral',verbose=False):
  #== set up outfile
  outfile = os.path.abspath(outfile)
  outdir = os.path.split(outfile)[0]
  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  if isinstance(segfile,str):
    segfile = segfile.split(',')
  if isinstance(segdeffile,str):
    segdeffile = segdeffile.split(',')

  #== customise plot
  pylab.rcParams.update({
                         "text.usetex": True,
                         #"font.family": 'serif',
                         "text.verticalalignment": "center",
                         #"lines.markersize": 12,
                         #"lines.markeredgewidth": 2,
                         "lines.linewidth": 5,
                         "font.size": 20,
                         "axes.titlesize": 20,
                         "axes.labelsize": 18,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16,
                         "legend.fontsize": 20})

  #== load triggers
  if verbose:
    print >>sys.stdout, "Loading triggers..."
  if isinstance(cache,str):
    cache = [cache]
  beforetrigs = []
  for f in cache:
    if f.endswith('xml'):
      beforetrigs += dqDataUtils.fromtrigxml(open(f,'r'),start,end,table)
    else:
      beforetrigs += dqDataUtils.fromtrigfile(open(f,'r'),start,end,etg)

  #== load segments
  if verbose:
    print >>sys.stdout, "Loading segments..."
  segs = segmentlist()
  if segfile:
    for f in segfile:
      f = os.path.realpath(f)
      if not os.path.isfile(f):	
        print >>sys.stdout, "Cannot find segment file "+f
        sys.exit(1)

      #== load segment definer file
      if segdeffile:
        idx = segfile.index(f)
        try:
          if segdeffile[idx].endswith('xml'):
            segdefs = dqSegmentUtils.fromsegmentxml(open(segdeffile[idx],'r'))
          else:
            segdefs = segmentsUtils.fromsegwizard(open(segdeffile[idx],'r'))
        except:
          segdefs = segmentlist([segment(start,end)])
      else:
        segdefs = segmentlist([segment(start,end)])

      if f.endswith('xml'):
        segs += dqSegmentUtils.fromsegmentxml(open(f,'r'))&segdefs
      else:
        segs += segmentsUtils.fromsegwizard(open(f,'r'))&segdefs


  #== set up vetoed/nonvetoed trigger lists
  if verbose:
    print >>sys.stdout, "Parsing vetoed triggers..."
  notvetoedtrigs = dqDataUtils.parse_trigs(beforetrigs,segs,inclusive=False)
  vetoedtrigs    = dqDataUtils.parse_trigs(beforetrigs,segs,inclusive=True)

  #== set plot time unit
  if verbose:
    print >>sys.stdout, "Plotting..."
  if (end-start) < 20000:
    time_axis_unit = 60
  elif (end-start) >= 20000 and (end-start) < 604800:
    time_axis_unit = 3600
  else:
    time_axis_unit = 86400
  time_unit = {60:'minutes',3600:'hours',86400:'days'}

  #== set up plot lists
  notvetoedtimes = [float(trig.get_peak()-start)/time_axis_unit\
                    for trig in notvetoedtrigs]
  vetoedtimes = [float(trig.get_peak()-start)/time_axis_unit\
                 for trig in vetoedtrigs]
  notvetoedsnr = [trig.snr for trig in notvetoedtrigs]
  vetoedsnr = [trig.snr for trig in vetoedtrigs]

  #== convert time
  

  #== generate plot: plot in legend markersize then rescale plot markers
  fig = pylab.figure(figsize=[12,6])
  ax  = fig.gca()
  m = 5
  p1 = ax.plot(notvetoedtimes,notvetoedsnr,'b.',label='Not vetoed',\
               markersize=3*m)
  p2 = ax.plot(vetoedtimes,vetoedsnr,'r.',label='Vetoed',\
               markersize=3*m)
  ax.semilogy()
  leg = ax.legend()
  for p in [p1,p2]:
    p = p[0]
    p.set_markersize(m)
  ax.set_ylim(5,max(notvetoedsnr+vetoedsnr)*1.01)
  ax.set_xlim(0,float(end-start)/time_axis_unit)
  start = datetime(*date.XLALGPSToUTC(start)[:6])\
              .strftime("%B %d %Y, %H:%M:%S %ZUTC")
  ax.set_xlabel('Time ('+time_unit[time_axis_unit]+') since '+str(start))
  ax.set_ylabel('SNR ('+etg+')')
  ax.set_title('Effect of '+vetoname+' on '+etg+' triggers')
  ax.grid(True,which='major')
  ax.grid(True,which='majorminor')
  fig.savefig(outfile)
  if verbose:
    print >>sys.stdout, "Done."

