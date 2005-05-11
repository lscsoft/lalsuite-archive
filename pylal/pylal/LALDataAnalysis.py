#!/usr/bin/env python
import sys, getopt
import matplotlib.cm
from matplotlib.patches     import Patch
from matplotlib.axes        import Axes
from matplotlib.collections import PolyCollection
from matplotlib.colors      import normalize, Colormap
from pylab    import *
from readMeta import *

class snglInspiral(readSnglInspiralTable,Axes): 

  def __init__(self, triggerfile):
    readSnglInspiralTable.__init__(self,triggerfile)

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
    gca().grid(True)

  def plot_snr_v_time(self):
    S3start = 751651213
    secsPerDay = 3600*24
    plot((self.end_time - S3start) / secsPerDay,self.snr,'rx')
    gca().grid(True)
 
  def histHistc_snr(self):
    subplot(211)  
    if len(self.snr) > 4:
      nbins = len(self.snr)/4
    else: nbins = 1
    nn,xx,patches = hist(self.snr,nbins,normed=0)
    gca().grid(True)
    nTrigsAbove = 0
    for i in range(0, len(nn)): 
      nTrigsAbove += nn[i]
    mm = []   
    for i in range(0, len(nn)):
      mm.append(nTrigsAbove)
      nTrigsAbove -= nn[i]

    for i in range (0, len(nn)):
      if nn[i] == 0: 
        nn[i] = 1
    lognn = log10( nn )
    logmm = log10( mm )

    subplot(212)
    bar(xx,logmm,(xx[-1]-xx[0])/nbins,color='r')
    bar(xx,lognn,(xx[-1]-xx[0])/nbins,color='g')
    gca().grid(True)
    #gca().set_yscale('log')
    #gca().set_ylim( (0.001,1000))


class doubleCoincInspiral(readSnglInspiralTable,Axes):

  def __init__(self, triggerfile):
    self.table1 = readSnglInspiralTable(triggerfile)
    self.table2 = readSnglInspiralTable(triggerfile)
    # can't do the following
    #readSnglInspiralTable.__init__(self.table1,triggerfile)


class SnglBurst(readSnglBurstTable,Axes,Patch,PolyCollection):

  def __init__(self, triggerfile):
    readSnglBurstTable.__init__(self,triggerfile)

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
                
    if shading == 'faceted': edgecolors =  (0,0,0,1), 
    else: edgecolors = 'None'
  
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
    if norm is not None: collection.set_clim(vmin, vmax)
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


def usage():
  print "LALDataAnalysis.py -x <xml-file>"
  sys.exit(0)

def main():
  # parse options
  try:
    opts,args=getopt.getopt(sys.argv[1:],"x:h",["xmlfile=","help"])            
  except getopt.GetoptError:
    usage()
    sys.exit(2)
  for o,a in opts:
    if o in ("-x","--xmlfile"):
             xmlfile = a
    if o in ("-h", "--help"):
             usage()
             sys.exit()
  
  # testing snglInspiral class
  trigs = snglInspiral(xmlfile)
  trigs.histHistc_snr()
  #trigs.plot_snr_v_time()
  #trigs.plot_snr_v_chisq()
  #trigs.summary()
  show()


# execute main if this module is explicitly invoked by the user
if __name__=="__main__":
        main()
