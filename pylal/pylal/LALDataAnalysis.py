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

  def nevents(self):
    return len(self.table)

  def mkarray(self, colname):
    myarray = asarray( [ self.table[i][colname] for i in range(self.nevents())] )
    return myarray

  def summary(self):
    subplot(221)
    snr = self.mkarray("snr")
    hist(snr)
    xlabel(r'SNR', size='x-large')
    ylabel(r'# triggers', size='x-large')
    subplot(222)
    mass1 = self.mkarray("mass1")
    hist(mass1)
    title(r'Excess power trigger')

  def plot_snr_v_chisq(self):
    plot(self.snr,self.chisq,'rx')
    title('SNR vs CHISQ')
    xlabel('snr')
    ylabel('chisq')
    gca().grid(True)

  def plot_snr_v_time(self):
    S3start = 751651213
    S4start = 793130413
    secsPerDay = 3600*24
    plot((self.end_time - S3start) / secsPerDay,self.snr,'rx')
    title('SNR vs TIME')
    xlabel('time')
    ylabel('snr')
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

  def __init__(self, triggerfile1, triggerfile2):
    self.table1 = readSnglInspiralTable(triggerfile1)
    self.table2 = readSnglInspiralTable(triggerfile2)
    # can't do the followingh
    #readSnglInspiralTable.__init__(self.table1,triggerfile)

  def plot_m1_v_m2(self):
    plot(self.table1.mass1,self.table1.mass2,'rx')
    plot(self.table2.mass1,self.table2.mass2,'b+')
    gca().grid(True)


class tripleCoincInspiral(readSnglInspiralTable,Axes):

  def __init__(self, triggerfile1, triggerfile2, triggerfile3):
    self.table1 = readSnglInspiralTable(triggerfile1)
    self.table2 = readSnglInspiralTable(triggerfile2)
    self.table3 = readSnglInspiralTable(triggerfile3)
    # can't do the followingh
    #readSnglInspiralTable.__init__(self.table1,triggerfile)

  def plot_m1_v_m2(self):
    plot(self.table1.mass1,self.table1.mass2,'r+')
    plot(self.table2.mass1,self.table2.mass2,'bx')
    plot(self.table3.mass1,self.table3.mass2,'gx')
    gca().grid(True)


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
  helpmsg = """
Usage: 
python LALDataAnalysis.py [OPTIONS] FILE1 FILE2 FILE3 ....

   -a, --snglInsp-snrVtime       plot snr vs time from single inspiral xml 
   -b, --snglInsp-snrVchisq      plot snr vs chisq from single inspiral xml
   -c, --snglInsp-histHistc-snr  plot snr histograms from single inspiral xml 
   -d, --snglInsp-summary        plot summary info from single inspiral xml
   -p, --show-plot               display plot 
   -s, --save-fig                save plot in .png and .ps format
   -t, --temporary-test		 only for developers to test this program
   -h, --help                    show this usage message

Example: 
python LALDataAnalysis.py -a -s G1-INSPVETO-SIRE.xml
...
Saved plot to file G1-INSPVETO-SIRE.snglInsp-snrVtime.png
Saved plot to file G1-INSPVETO-SIRE.snglInsp-snrVtime.ps
  """
  print helpmsg


def main():
  # initialise variables
  show_plot = False
  save_fig  = False
  xml_file  = None
  plot_type = None

  # define short and long option names
  shortopts = "abcdpsth"
  longopts = [
    "snglInsp-snrVtime",
    "snglInsp-snrVchisq",
    "snglInsp-histHistc-snr",
    "snglInsp-summary",
    "show-plot",
    "save-fig",
    "temporary-test",
    "help"]
  
  # parse command line options
  try:
    options,arguments = getopt.getopt(sys.argv[1:],shortopts,longopts)            
  except getopt.GetoptError:
    usage()
    sys.exit(2)

  print options
  print arguments

  for o,a in options:
    if o in ("-a","--snglInsp-snrVtime"):
      xml_files  = arguments
      plot_type = "snglInsp-snrVtime"
    if o in ("-b","--snglInsp-snrVchisq"):
      xml_files  = arguments
      plot_type = "snglInsp-snrVchisq"
    if o in ("-c","--snglInsp-histHistc-snr"):
      xml_files  = arguments
      plot_type = "snglInsp-histHistc-snr"
    if o in ("-d","--snglInsp-summary"):
      xml_files  = arguments
      plot_type = "snglInsp-summary"
    if o in ("-p","--show-plot"):
      show_plot = True 
    if o in ("-s","--save-fig"):
      save_fig  = True   
    if o in ("-t","--temporary-test"):
      # test new code here
      sys.exit(0)
    if o in ("-h", "--help"):
      usage()
      sys.exit(0)

  # check options are sane 
  if not plot_type:
    print >> sys.stderr, "No plot option specified"
    sys.exit(1)
  if not xml_files:
    print >> sys.stderr, "No trigger file specified"
    sys.exit(1)

  # read data file and call plotting function desired
  if plot_type == "snglInsp-snrVtime":
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_time()
  elif plot_type == "snglInsp-snrVchisq":
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_chisq()
  elif plot_type == "snglInsp-histHistc-snr":
    trigs = snglInspiral(xml_files[0])
    trigs.histHistc_snr()
  elif plot_type == "snglInsp-summary":
    trigs = snglInspiral(xml_files[0])
    trigs.summary()
  

  # save and show plot if desired
  if save_fig:
    png_file = xml_file[:-3] + plot_type + ".png"
    ps_file  = xml_file[:-3] + plot_type + ".ps"
    savefig(png_file)
    savefig(ps_file)
    print "Saved plot to file %s" % (png_file)
    print "Saved plot to file %s" % (ps_file)
  if show_plot:
    show()

# execute main if this module is explicitly invoked by the user
if __name__=="__main__":
        main()
