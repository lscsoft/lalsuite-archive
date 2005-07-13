#!/usr/bin/env python
import sys, getopt
import matplotlib.cm
from matplotlib.patches     import Patch
from matplotlib.axes        import Axes
from matplotlib.collections import PolyCollection
from matplotlib.colors      import normalize, Colormap
from optparse import * 
from pylab    import *
from readMeta import *

def simpleplot(*args):
  if len(args)==3:
      mytable, col1, col2 = args
  else:
      raise TypeError, 'Illegal arguments to simpleplot; see help(rectfill)'
   
  if mytable is not None: assert(isinstance(mytable, metaDataTable))
  tmpvar1 = mytable.mkarray(col1)
  tmpvar2 = mytable.mkarray(col2)
  plot(tmpvar1,tmpvar2,'rx')
  xlabel(col1, size='x-large')
  ylabel(col2, size='x-large')
  title(col1 + ' v ' + col2)
  grid(True)



def tfplot(*args, **kwargs):
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


def main():
  # define usage and command line options and arguments - parse
  usage = "usage: %prog ..."
  parser = OptionParser( usage )

  opts_snglInsp = OptionGroup( parser, "Single Inspiral plotting functions",\
	"Example ..." )
  opts_snglInsp.add_option( "-a", "--snglInsp_snrVtime",\
	action="store_true", default=False,\
	help="plot snr vs time from a single inspiral xml" )
  opts_snglInsp.add_option( "-b", "--snglInsp_snrVchisq",\
        action="store_true", default=False,\
        help="plot snr vs chisq from single inspiral xml")
  opts_snglInsp.add_option( "-c", "--snglInsp_histHistc_snr",\
        action="store_true", default=False,\
        help="plot snr histograms from single inspiral xml" )
  opts_snglInsp.add_option( "-d", "--snglInsp_summary",\
        action="store_true", default=False,\
        help="plot summary info from single inspiral xml" )

  parser.add_option_group( opts_snglInsp )

  parser.add_option( "-p", "--show_plot",\
        action="store_true", default=False,\
        help="display plot" )
  # change this so that by default the fig is always saved using 
  # the name convention already implemented. Now instead of --save-fig
  # you have --save-off and --save-as flag to override
  # the standard name. Also add the date+time to the standard name OR
  # check for the existence of the standard name + 001, 002, 003, ...
  # Remove where possible the defns of dest in favour of the long option name
  parser.add_option( "-s", "--save_fig",\
        action="store_true", default=False,\
        help="save figure in .png and .ps format" )
  parser.add_option( "-t", "--temporary-test",\
        action="store_true", default=False,\
        help="only for developers to test this program" )

  ( options, xml_files ) = parser.parse_args()
  
  # check xml_files have been given as required arguments 
  if not xml_files:
    print >> sys.stderr, "No trigger file specified"
    sys.exit(1)

  # read data files and call plotting function desired
  if   options.snglInsp_snrVtime:
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_time()
  elif options.snglInsp_snrVchisq:
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_chisq()
  elif options.snglInsp_histHistc_snr:
    trigs = snglInspiral(xml_files[0])
    trigs.histHistc_snr()
  elif options.snglInsp_summary:
    trigs = snglInspiral(xml_files[0])
    trigs.summary()
  else:
    print >> sys.stderr, "No plot option specified"
    sys.exit(1)
  
  # save and show plot if desired
  if options.save_fig:
    png_file = xml_file[:-3] + plot_type + ".png"
    ps_file  = xml_file[:-3] + plot_type + ".ps"
    savefig(png_file)
    savefig(ps_file)
    print "Saved plot to file %s" % (png_file)
    print "Saved plot to file %s" % (ps_file)
  if options.show_plot:
    show()

# execute main if this module is explicitly invoked by the user
if __name__=="__main__":
        main()
