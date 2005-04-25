#!/usr/bin/env python
import sys
from matplotlib.patches import Patch
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.colors import normalize, Colormap
import matplotlib.cm
from pylab import *
from lgen import metaio

class SnglBurstTable(Axes,Patch,PolyCollection): 
  
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sngl_burst('%s' % triggerfile)
    self.start_time = M[:,0]+1.0e-9*M[:,1]
    self.peak_time = M[:,2]+1.0e-9*M[:,3]
    self.duration = M[:,4]
    self.central_freq = M[:,5]
    self.bandwidth = M[:,6]
    self.snr = M[:,7]
    self.confidence = M[:,8]
    #subplot(221)
    #self.tfplane=gca()
    

  def compute_tf_volume(self):
    return self.duration * self.bandwidth

  def histogram_confidence(self):
    hist(self.confidence)

  def tfplot(self,  *args, **kwargs):
        """\
tfplot(x, y, s=20, c='b', marker='o', cmap=None, norm=None,
        vmin=None, vmax=None, alpha=1.0)

Supported function signatures:

    TFPLOT(x, y) - make a scatter plot of x vs y

    TFPLOT(x, y, s) - make a scatter plot of x vs y with size in area
      given by s

    TFPLOT(x, y, s, c) - make a scatter plot of x vs y with size in area
      given by s and colors given by c

    TFPLOT(x, y, s, c, **kwargs) - control colormapping and scaling
      with keyword args; see below

Make a scatter plot of x versus y.  s is a size in points^2 a scalar
or an array of the same length as x or y.  c is a color and can be a
    """
        shading = kwargs.get('shading', 'faceted')
        cmap = kwargs.get('cmap', cm.get_cmap())
        norm = kwargs.get('norm', normalize())
        alpha = kwargs.get('alpha', 1.0)
        vmin = kwargs.get('vmin', None)
        vmax = kwargs.get('vmax', None)  

        if len(args)==5:
            X, dX, Y, dY, C = args
        else:
            raise TypeError, 'Illegal arguments to rectfill; see help(rectfill)'
       
        Nx, = X.shape

        verts = [ ( (X[i,] , Y[i,]) , (X[i,]+dX[i,] , Y[i,]),
                    (X[i,]+dX[i,] , Y[i,]+dY[i,]), 
                    (X[i,] , Y[i,]+dY[i,]) )
                  for i in range(Nx-1) ]
           
        C = array([C[i,] for i in range(Nx-1)])

                
        if shading == 'faceted':
            edgecolors =  (0,0,0,1), 
        else:
            edgecolors = 'None'

        collection = PolyCollection(
            verts,
            edgecolors   = edgecolors,
            antialiaseds = (0,),
            linewidths   = (0.25,),
            )

        collection.set_alpha(alpha)
        collection.set_array(C)
        if norm is not None: assert(isinstance(norm, normalize))
        if cmap is not None: assert(isinstance(cmap, Colormap))
        collection.set_cmap(cmap)
        collection.set_norm(norm)

        if norm is not None:
            collection.set_clim(vmin, vmax)

        a = gca()
        a.grid(False)

        minx = amin(X)
        maxx = amax(X)
        miny = amin(Y)
        maxy = amax(Y)

        corners = (minx, miny), (maxx, maxy) 
        
        a.update_datalim( corners )
        a.autoscale_view()

        # add the collection last
        a.add_collection(collection)
        xlabel(r'Time (secs)')
        ylabel(r'Frequency (Hz)')

        return collection


class SnglInspiralTable(Axes): 
  
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sngl_inspiral('%s' % triggerfile)
    self.end_time = M[:,0]+1.0e-9*M[:,1]
    self.eff_distance = M[:,2]
    self.coa_phase = M[:,3]
    self.mass1 = M[:,4]
    self.mass2 = M[:,5]
    self.mchirp = M[:,6]
    self.eta = M[:,7]
    self.snr = M[:,8]
    self.chisq = M[:,9]

  def summary(self):
    subplot(221)
    hist(self.snr)
    xlabel(r'SNR', size='x-large')
    ylabel(r'# triggers', size='x-large')
    subplot(222)
    hist(self.mass1)
    title(r'Excess power trigger')

  def plot_snr_v_chisq(self):
    plot(self.snr,self.chisq)


class SimInspiralTable(Axes): 
  
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sim_inspiral('%s' % triggerfile)
    self.geo_end_time = M[:,0]+1.0e-9*M[:,1]
    self.distance = M[:,2]
    self.f_final = M[:,3]
    self.mass1 = M[:,4]
    self.mass2 = M[:,5]
    self.mchirp = M[:,6]
    self.eta = M[:,7]
    self.eff_dist_h = M[:,8]
    self.eff_dist_l = M[:,9]

  def summary(self):
    subplot(221)
    hist(self.distance)
    xlabel(r'Distance', size='x-large')
    ylabel(r'# injections', size='x-large')
    subplot(222)
    hist(self.eff_dist_h,20)
    title(r'Excess power trigger')
