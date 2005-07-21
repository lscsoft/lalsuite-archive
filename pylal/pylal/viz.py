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

##############################################
# function to read in a column from two tables 
def readcolfrom2tables(table1, table2, col_name ):
  
  if table1 is not None: assert(isinstance(table1, metaDataTable))
  if table2 is not None: assert(isinstance(table2, metaDataTable))

  if table1.nevents() != table2.nevents():
    print >>sys.stderr, "nevents in table1 and table2 must be equal"
    sys.exit(1)
 
  if table1.table[1].has_key('ifo'):
    ifo = table1.table[1]["ifo"]
  elif table2.table[1].has_key('ifo'):
    ifo = table2.table[1]["ifo"]
  else:
    ifo = None

  col_names = [col_name]
  if ifo:
    col_names.append(col_name + '_' + ifo[0].lower())
    col_names.append(ifo[0].lower() + '_' + col_name)

  if 'dist' in col_name:
    col_names.append(col_name + 'ance')

  for c_name in col_names:
    
    if table1.table[1].has_key(c_name):
      tmpvar1 = table1.mkarray(c_name)
      if 'time' in c_name:
        tmpvar1 = tmpvar1 + 1e-9 * table1.mkarray(c_name + '_ns')

    if table2.table[1].has_key(c_name):
      tmpvar2 = table2.mkarray(c_name)
      if 'time' in c_name:
        tmpvar2 = tmpvar2 + 1e-9 * table2.mkarray(c_name + '_ns')

  cols = []
  cols.append(tmpvar1)
  cols.append(tmpvar2)
  cols.append(ifo)

  return cols

  
#################################################
# function to read in a column from a given table 
def readcol(table, col_name, ifo=None ):
  
  if table is not None: assert(isinstance(table, metaDataTable))

  col_names = [col_name]
  if ifo:
    col_names.append(col_name + '_' + ifo[0].lower())
    col_names.append(ifo[0].lower() + '_' + col_name)

  if 'dist' in col_name:
    col_names.append(col_name + 'ance')

  for c_name in col_names:
    
    if table.table[1].has_key(c_name):
      col_data = table.mkarray(c_name)
      if 'time' in c_name:
        col_data = col_data + 1e-9 * table.mkarray(c_name + '_ns')


  return col_data
 
###################################################
# function to plot the difference between values of 'col_name' in
# two tables, table1 and table2
def plotdiff(table1, table2, col_name, plot_type, plot_sym, units=None, 
  axis_range=[0,0,0,0], output_name = None):
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  tmp_diff = tmpvar2 - tmpvar1

  if plot_type == 'plot':
    plot(tmpvar1, tmp_diff, plot_sym, markersize=12)
  elif plot_type == 'log':
    semilogx(tmpvar1, tmp_diff, plot_sym, markersize=12)
    
  if units:
    xlabel(col_name + ' (' + units +')', size='x-large')
    ylabel(col_name + ' difference (' + units +')', size='x-large')
  else:
    xlabel(col_name, size='x-large')
    ylabel(col_name + ' difference', size='x-large')
  
  if not axis_range[0]:
    axis_range[0] = min(tmpvar1)
  if not axis_range[1]:
    axis_range[1] = max(tmpvar1)
  if not axis_range[2]:
    if min(tmp_diff) > 0:
      axis_range[2] = 0.9 * min(tmp_diff)
    else:
      axis_range[2] = 1.1 * min(tmp_diff)
  if not axis_range[3]:
    if max(tmp_diff) > 0:
      axis_range[3] = 1.1 * max(tmp_diff)
    else:
      axis_range[3] = 0.9 * max(tmp_diff)

  if ifo:
    title(ifo + ' ' + col_name + '  Accuracy', size='x-large')
  else:
    title(col_name + '  Accuracy', size='x-large')
  
  
  axis(axis_range)
  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name' in
# two tables, table1 and table2
def plotfracdiff(table1, table2, col_name, plot_type, plot_sym, units=None, 
  axis_range=[0,0,0,0], output_name = None):
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  frac_diff = (tmpvar2 - tmpvar1)/tmpvar1

  if plot_type == 'plot':
    plot(tmpvar1, frac_diff,plot_sym,markersize=12)
  elif plot_type == 'log':
    semilogx(tmpvar1, frac_diff,plot_sym,markersize=12)
    
  if units:
    xlabel(col_name + ' (' + units +')', size='x-large')
  else:
    xlabel(col_name, size='x-large')
  
  ylabel(col_name + ' fractional difference', size='x-large')
  
  if not axis_range[0]:
    axis_range[0] = min(tmpvar1)
  if not axis_range[1]:
    axis_range[1] = max(tmpvar1)
  if not axis_range[2]:
    if min(frac_diff) > 0:
      axis_range[2] = 0.9 * min(frac_diff)
    else:
      axis_range[2] = 1.1 * min(frac_diff)
  if not axis_range[3]:
    if max(frac_diff) > 0:
      axis_range[3] = 1.1 * max(frac_diff)
    else:
      axis_range[3] = 0.9 * max(frac_diff)

  if ifo:
    title(ifo + ' ' + col_name + '  Accuracy', size='x-large')
  else:
    title(col_name + '  Accuracy', size='x-large')
  
  
  axis(axis_range)
  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_frac_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name_a' in
# two tables, table1 and table2 against the values of 'col_name_b' in table1
def plotdiffa_vs_b(table1, table2, col_name_a, col_name_b, plot_type, 
  plot_sym, units_a=None, units_b=None,axis_range=[0,0,0,0], 
  output_name = None):
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name_a)

  diff_a = (tmpvar2 - tmpvar1)
  col_b = readcol(table1, col_name_b, ifo ) 

  if plot_type == 'plot':
    plot(col_b, diff_a,plot_sym,markersize=12)
  elif plot_type == 'log':
    semilogx(col_b, diff_a,plot_sym,markersize=12)
    
  if units_b:
    xlabel(col_name_b + ' (' + units_b +')', size='x-large')
  else:
    xlabel(col_name_b, size='x-large')
  
  if units_a:
    ylabel(col_name_a + ' difference (' + units_a + ')', size='x-large')
  else:
    ylabel(col_name_a + ' difference', size='x-large')
    
  if not axis_range[0]:
    axis_range[0] = min(col_b)
  if not axis_range[1]:
    axis_range[1] = max(col_b)
  if not axis_range[2]:
    if min(diff_a) > 0:
      axis_range[2] = 0.9 * min(diff_a)
    else:
      axis_range[2] = 1.1 * min(diff_a)
  if not axis_range[3]:
    if max(diff_a) > 0:
      axis_range[3] = 1.1 * max(diff_a)
    else:
      axis_range[3] = 0.9 * max(diff_a)

  if ifo:
    title(ifo + ' ' + col_name_a + '  Accuracy', size='x-large')
  else:
    title(col_name + '  Accuracy', size='x-large')
  
  
  axis(axis_range)
  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '_accuracy.png'
    savefig(output_name)
 



######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histdiff(table1, table2, col_name, plot_type, sym, units=None, 
  nbins = None, width = None, output_name = None):
  
  histcolors = ['r','b','k']
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  tmp_diff = tmpvar2 - tmpvar1

  if (plot_type == 'frac_hist'):
    tmp_diff /= tmpvar1

  if not nbins:
    nbins = 10
  
  bins = []
  if width:
    for i in range(-nbins,nbins):
      bins.append(width * i/nbins)
  
  if bins:
    out = hist(tmp_diff,bins)
  else:
    out = hist(tmp_diff,nbins)
 
  width = out[1][1] - out[1][0]
  bar(out[1],out[0],width,color=histcolors[sym])

  figtext(0.13,0.8 - 0.1* sym," mean = %6.3e" % mean(tmp_diff))
  figtext(0.13,0.75 - 0.1 * sym,'sigma = %6.3e' % std(tmp_diff))
   
  label = col_name 
  if (plot_type == 'frac_hist'):
    label += ' fractional'
  label += ' difference'
  if units and not (plot_type == 'frac_hist'):
    label += ' (' + units +')'
  xlabel(label, size='x-large')

  ylabel('Number', size='x-large')
  
  if ifo:
    title(ifo + ' ' + col_name + '  Histogram', size='x-large')
  else:
    title(col_name + ' Histogram', size='x-large')
  
  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo  
    if (plot_type == 'frac_hist'):  
      output_name += '_' + col_name + '_frac_histogram.png'
    else:
      output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)
  

###################################################
# function to plot the value of 'col_name' in table1 vs its
# value in table2 
def plotval(table1, table2, col_name, plot_type, units=None, xlab=None, \
  ylab=None, axis_range=[0,0,0,0], output_name=None):
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  if plot_type == 'plot':
    plot(tmpvar1, tmpvar2,'bx')
  elif plot_type == 'log':
    loglog(tmpvar1, tmpvar2,'bx')
    
  if xlab:
    xlab += ' ' + col_name
  else:
    xlab = col_name

  if ylab:
    ylab += ' ' + col_name
  else:
    ylab = col_name
    
  if units:
    xlab +=' (' + units +')'
    ylab += ' (' + units +')'
    
  xlabel(xlab, size='x-large')
  ylabel(ylab, size='x-large')
 
  
  if not axis_range[0]:
    axis_range[0] = min(tmpvar1)
  if not axis_range[1]:
    axis_range[1] = max(tmpvar1)

  if not axis_range[0]:
    if min(tmpvar1) > 0:
      axis_range[0] = 0.9 * min(tmpvar1)
    else:
      axis_range[0] = 1.1 * min(tmpvar1)
  if not axis_range[1]:
    if max(tmpvar1) > 0:
      axis_range[1] = 1.1 * max(tmpvar1)
    else:
      axis_range[1] = 0.9 * max(tmpvar1)
  
  if not axis_range[2]:
    if min(tmpvar2) > 0:
      axis_range[2] = 0.9 * min(tmpvar2)
    else:
      axis_range[2] = 1.1 * min(tmpvar2)
  if not axis_range[3]:
    if max(tmpvar2) > 0:
      axis_range[3] = 1.1 * max(tmpvar2)
    else:
      axis_range[3] = 0.9 * max(tmpvar2)

  if ifo:
    title(ifo + ' ' + col_name, size='x-large')
  else:
    title(col_name, size='x-large')
  
  
  axis(axis_range)
  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_plot.png'
    savefig(output_name)
 

###################################################
# function to plot the difference between values of 'col_name' in
# two tables, table1 and table2
def histdiffdiff(ifo1_trig, ifo2_trig, inj, col_name, sym, units=None, 
  nbins = None, width = None, output_name = None):
 
  histcolors = ['b','r','k']

  [tmpvar1, injvar1, ifo1 ] = readcolfrom2tables(ifo1_trig, inj, col_name)
  [tmpvar2, injvar2, ifo2 ] = readcolfrom2tables(ifo2_trig, inj, col_name)
  
  diff1 = injvar1 - tmpvar1
  diff2 = injvar2 - tmpvar2
  
  diffdiff = diff1 - diff2
  
  if not nbins:
    nbins = 10
  
  bins = []
  if width:
    for i in range(-nbins,nbins):
      bins.append(width * i/nbins)

  if bins:
    out = hist(diffdiff,bins)
  else:
    out = hist(diffdiff,nbins)

  width = out[1][1] - out[1][0]
  bar(out[1],out[0],width,color=histcolors[sym])

  figtext(0.13,0.8 - 0.1* sym," mean = %6.3e" % mean(diffdiff))
  figtext(0.13,0.75 - 0.1 * sym,'sigma = %6.3e' % std(diffdiff))
 
  label = col_name 
  label += ' difference'
  if units:
    label += ' (' + units +')'
  xlabel(label, size='x-large')

  ylabel('Number', size='x-large')
  
  title(ifo1 + ' - ' + ifo2 + ' ' + col_name + ' difference histogram', \
    size='x-large')
  
  grid(True)

  if output_name:
    output_name += '_' + ifo1 + '_' + ifo2 
    output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)
  

   

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
