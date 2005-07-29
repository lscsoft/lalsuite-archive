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


#################################################
# function to read in a column from a given table 
def readcol(table, col_name, ifo=None ):
  """
  function to read in a column from a given table.  If the column is ifo
  dependent (eg, end_time or eff_dist) then the ifo is used to select the
  appropriate value.
  
  @param table: metaDataTable
  @param col_name: name of column to read in 
  @param ifo: name of ifo (used to extract the appropriate end_time/eff_dist
  """
  
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

##############################################
# function to read in a column from two tables 
def readcolfrom2tables(table1, table2, col_name ):
  """
  function to read in a column from a two tables.  If the column is ifo
  dependent (eg, end_time or eff_dist) then the ifo of one table is used to 
  select the appropriate value in the other table.  The number of events in
  the two tables must be equal.
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to read in 
  """
  
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

  col1 = readcol(table1, col_name, ifo)
  col2 = readcol(table2, col_name, ifo)

  cols = []
  cols.append(col1)
  cols.append(col2)
  cols.append(ifo)

  return cols


##############################################
# function to read in a column from two tables 
def timeindays(col_data ):
  """
  function to re-express the time in days after the start of the run
  
  @param col_data: array containing times
  """
  s2times = [729273613, 734367613]
  s3times = [751658413, 757699213]
  s4times = [793130413, 795679213]
  
  if col_data[0] > s2times[0] and col_data[0] < s2times[1]:
    start = s2times[0]
  elif col_data[0] > s3times[0] and col_data[0] < s3times[1]:
    start = s3times[0]
  elif col_data[0] > s4times[0] and col_data[0] < s4times[1]:
    start = s4times[0]
  else:
    print >> sys.stderr, "events not from a known science run"

  col_data = (col_data - start)/(60 * 60 * 24.0)

  return col_data

#################################################################
# function to plot the col1 vs col2 from the table
def plot_a_v_b(table, col_name_a, col_name_b, plot_type, plot_sym, \
  output_name = None):

  if table.table[1].has_key('ifo'):
    ifo = table.table[1]["ifo"]
  else:
    ifo = None

  col_a = readcol(table, col_name_a, ifo )
  col_b = readcol(table, col_name_b, ifo )

  if 'time' in col_name_a:
    col_a = timeindays(col_a)
  if 'time' in col_name_b:
    col_b = timeindays(col_b)
    
  if plot_type == 'linear':
    plot(col_a, col_b,plot_sym, markersize=12,markerfacecolor=None)
  elif plot_type == 'logx':
    semilogx(col_a, col_b, plot_sym, markersize=12,markerfacecolor=None)
  elif plot_type == 'logy':
    semilogx(col_a, col_b, plot_sym, markersize=12,markerfacecolor=None)
  elif plot_type == 'loglog':
    loglog(col_a, col_b, plot_sym, markersize=12,markerfacecolor=None)
  
  if ifo:
    title(ifo + ' ' + col_name_a + ' vs ' + col_name_b, size='x-large')
  else:
    title(col_name_a + ' vs ' + col_name_b, size='x-large')

  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '.png'
    savefig(output_name)

#################################################################
# function to plot the difference between values of 'col_name' in
# two tables, table1 and table2
def plotdiff(table1, table2, col_name, plot_type, plot_sym):
  """
  function to plot the difference between the value of col_name stored in 2
  tables (of equal length).  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  tmp_diff = tmpvar2 - tmpvar1

  if 'time' in col_name:
    tmpvar1 = timeindays(tmpvar1)

  if plot_type == 'linear':
    plot(tmpvar1, tmp_diff, plot_sym, markersize=12,markerfacecolor=None)
  elif plot_type == 'log':
    semilogx(tmpvar1, tmp_diff, plot_sym, markersize=12,markerfacecolor=None)
    
#################################################################
# function to label above plot
def labeldiff(col_name, units = None, axis = [0,0,0,0], leg = None, 
  title_text = None, output_name = None):
  """
  function to label the output of plotdiff
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """
  
  if units:
    xlabel(col_name + ' (' + units +')', size='x-large')
    ylabel(col_name + ' difference (' + units +')', size='x-large')
  else:
    xlabel(col_name, size='x-large')
    ylabel(col_name + ' difference', size='x-large')

  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name + '  Accuracy', size='x-large',
      weight='bold')
  else:
    title(col_name + ' Accuracy', size='x-large',weight='bold')
  
  if output_name:
    output_name += '_' + col_name + '_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name' in
# two tables, table1 and table2
def plotfracdiff(table1, table2, col_name, plot_type, plot_sym):
  """
  function to plot the fractional difference between the value of 
  col_name stored in 2 tables (of equal length).  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  frac_diff = (tmpvar2 - tmpvar1)/tmpvar1

  if 'time' in col_name:
    tmpvar1 = timeindays(tmpvar1)

  if plot_type == 'linear':
    plot(tmpvar1, frac_diff,plot_sym,markersize=12,markerfacecolor=None)
  elif plot_type == 'log':
    semilogx(tmpvar1, frac_diff,plot_sym,markersize=12,markerfacecolor=None)


#################################################################
# function to label above plot
def labelfracdiff(col_name, units = None, axis = [0,0,0,0], leg = None, 
  title_text = None, output_name = None):
  """
  function to label the output of plotfracdiff
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """

  if units:
    xlabel(col_name + ' (' + units +')', size='x-large')
  else:
    xlabel(col_name, size='x-large')
  
  ylabel(col_name + ' fractional difference', size='x-large')

  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name + '  Accuracy', size='x-large',
      weight='bold')
  else:
    title(col_name + ' Accuracy', size='x-large',weight='bold')
  
  if output_name:
    output_name += '_' + col_name + '_frac_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name_a' in
# two tables, table1 and table2 against the values of 'col_name_b' in table1
def plotdiffa_vs_b(table1, table2, col_name_a, col_name_b, plot_type, 
  plot_sym):
  """
  function to plot the difference if col_name_a in two tables against the
  value of col_name_b in table1.  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name_a: name of column to plot difference of on y-axis
  @param col_name_b: name of column to plot on x-axis
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name_a)

  diff_a = (tmpvar2 - tmpvar1)
  col_b = readcol(table1, col_name_b, ifo ) 

  if 'time' in col_name_b:
    col_b = timeindays(col_b)

  if plot_type == 'linear':
    plot(col_b, diff_a,plot_sym,markersize=12,markerfacecolor=None)
  elif plot_type == 'log':
    semilogx(col_b, diff_a,plot_sym,markersize=12,markerfacecolor=None)
  
  
#################################################################
# function to label above plot
def labeldiffa_vs_b(col_name_a, col_name_b, units_a=None, units_b=None,
  axis=[0,0,0,0], leg = None, title_text = None, output_name = None):
  """
  function to label the output of plotdiffa_vs_b
  
  @param col_name_a: name of column to plot
  @param col_name_b: name of column to plot
  @param units_a: the units of column a
  @param units_b: the units of column b
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """
    
  if units_b:
    xlabel(col_name_b + ' (' + units_b +')', size='x-large')
  else:
    xlabel(col_name_b, size='x-large')
  
  if units_a:
    ylabel(col_name_a + ' difference (' + units_a + ')', size='x-large')
  else:
    ylabel(col_name_a + ' difference', size='x-large')
 
  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name_a + ' Accuracy vs ' + col_name_b, 
      size='x-large', weight='bold')
  else:
    title(col_name_a + ' Accuracy vs ' + col_name_b,
      size='x-large',weight='bold')

  if output_name:
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '_accuracy.png'
    savefig(output_name)
 
###################################################
# function to plot the value of 'col_name' in table1 vs its
# value in table2 
def plotval(table1, table2, col_name, plot_type, plot_sym):
  """
  function to plot the value of col_name in table1 vs its value in table2  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  if plot_type == 'linear':
    plot(tmpvar1, tmpvar2,plot_sym,markersize=12,markerfacecolor=None)
  elif plot_type == 'log':
    loglog(tmpvar1, tmpvar2,plot_sym,markersize=12,markerfacecolor=None)
   

#################################################################
# function to label above plot
def labelval(col_name, units = None, axis = [0,0,0,0], xlab = None, \
  ylab = None, leg = None, title_text = None, output_name = None):
  """
  function to label the output of plotval
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param xlab: label for x-axis
  @param ylab: label for y-axis
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """

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
 
  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name, size='x-large', weight='bold')
  else:
    title(col_name,size='x-large',weight='bold')

  if output_name:
    output_name += '_' + col_name + '_plot.png'
    savefig(output_name)
  

###########################################################
# function to plot the value of 'col_name' in ifo1  vs its
# value in ifo2 
def plotcoincval(coinctable, col_name, ifo1, ifo2, plot_sym, plot_type):
  
  ifo_coinc = coinctable.coinctype(ifo1,ifo2)

  ifo1_val = ifo_coinc.mkarray(col_name,ifo1)
  ifo2_val = ifo_coinc.mkarray(col_name,ifo2)
  
  if plot_type == 'linear':
    plot(ifo1_val, ifo2_val,plot_sym,markersize=12)
  elif plot_type == 'log':
    loglog(ifo1_val, ifo2_val,plot_sym,markersize=12)

###########################################################
# function to plot the value of 'col_name' in hanford  vs its
# value in ifo.  
def plotcoinchanford(coinctable, col_name, ifo, hanford_method, plot_sym,\
  plot_type):
  
  ifo_coinc = coinctable.coinctype(ifo,'H1','H2')

  ifo_val = ifo_coinc.mkarray(col_name,ifo)
  h1_val = ifo_coinc.mkarray(col_name,'H1')
  h2_val = ifo_coinc.mkarray(col_name,'H2')

  if hanford_method == 'sum':
    h_val = h1_val + h2_val
  if hanford_method == 'mean':
    h_val = (h1_val + h2_val)/2
   
  if ifo_val: 
    if plot_type == 'linear':
      plot(ifo_val, h_val,plot_sym,markersize=12)
    elif plot_type == 'log':
      loglog(ifo1_val, ifo2_val,plot_sym,markersize=12)




######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histcol(table1, col_name,nbins = None, width = None, output_name = None):
 
  if table1.table[1].has_key('ifo'):
    ifo = table1.table[1]["ifo"]
  else:
    ifo = None

  data = readcol(table1, col_name, ifo )

  if not nbins:
    nbins = 10
  
  bins = []
  if width:
    for i in range(-nbins,nbins):
      bins.append(width * i/nbins)
  
  xlabel(col_name, size='x-large')
  ylabel('Number', size='x-large')

  if bins:
    out = hist(data,bins)
  else:
    out = hist(data,nbins)

  width = out[1][1] - out[1][0]
  bar(out[1],out[0],width)
  
  if ifo:
    title(ifo + ' ' + col_name + ' histogram', size='x-large')
  else:
    title(col_name + ' histogram', size='x-large')

  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)


######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def cumhistcol(table1, col_name, normalization=None,output_name = None):
 
  if table1.table[1].has_key('ifo'):
    ifo = table1.table[1]["ifo"]
  else:
    ifo = None

  data = readcol(table1, col_name, ifo )

  data_len = len(data)
  data_sort = sort(data)
  data_range = arange(len(data))

  y_data = data_len - data_range
  if normalization:
    y_data = y_data/float(normalization)
  
  semilogy(data_sort, y_data,'b-')
  
  xlabel(col_name, size='x-large')
  
  if normalization:
    ylabel('Probability', size='x-large')
  else:  
    ylabel('Cumulative Number', size='x-large')

  if ifo:
    title_string = ifo + ' ' + col_name
  else:
    title_string = col_name
  if normalization:
    title_string += ' normalized'
  title_string += ' cumulative histogram'
  title(title_string, size='x-large')

  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    if normalization:
      output_name += '_' + col_name + '_norm_hist.png'
    else:
      output_name += '_' + col_name + '_cum_hist.png'
    savefig(output_name)


######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def cumhistsnr(trigs, ifos = None, min_val = None, max_val = None, \
  nbins = None):
 
  if ifos:
    trigs = trigs.coinctype(ifos[0],ifos[1])
  snr = asarray([ pow(trigs.table[i]["snrsq"],0.5) for i in range(trigs.nevents()) ] )

  # set up the bin boundaries
  if not nbins:
    nbins = 20


  bins = []
  if max_val and min_val:
    for i in range(nbins):
      bins.append(min_val + i*(max_val - min_val)/nbins)
  
  if bins:
    [num,bin,info] = hist(snr,bins)
  else:
    [num,bin,info] = hist(snr,nbins)
 
  cum_num = [sum(num)]
  cum_num.extend(sum(num) - cumsum(num))
  cum_num.pop()
  figure(23)
  semilogy(cum_num,bins,'kx')



######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histdiff(table1, table2, col_name, plot_type, hist_col, nbins=None, 
  width=None):
  """
  function to plot a histogram of the difference of the value of col_name
  between table1 and table2  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'hist' or 'frac_hist' 
  @param hist_col: the colour of the histogram
  @param nbins: number of bins to plot in histogram (default = 10)
  @param width: the maximum difference to be shown
  """

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
  bar(out[1],out[0],width,color=hist_col)

  #figtext(0.13,0.8 - 0.1* sym," mean = %6.3e" % mean(tmp_diff))
  #figtext(0.13,0.75 - 0.1 * sym,'sigma = %6.3e' % std(tmp_diff))
  
######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def labelhistdiff(col_name, plot_type, units, leg = None, title_text = None, 
  output_name = None):
  """
  function to label the output of histdiff

  @param col_name: name of column to plot
  @param plot_type: either 'hist' or 'frac_hist' 
  @param units: the units of the column
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """


  label = col_name 
  if (plot_type == 'frac_hist'):
    label += ' fractional'
  label += ' difference'
  if units and not (plot_type == 'frac_hist'):
    label += ' (' + units +')'
  xlabel(label, size='x-large')

  ylabel('Number', size='x-large')
  
  if title_text:
    title(title_text + ' ' + col_name + '  Histogram', size='x-large')
  else:
    title(col_name + ' Histogram', size='x-large')
  
  grid(True)

  if output_name:
    if (plot_type == 'frac_hist'):  
      output_name += '_' + col_name + '_frac_histogram.png'
    else:
      output_name += '_' + col_name + '_histogram.png'
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
  

######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histslides(slide_trigs, zerolag_trigs = None, ifos = None):
  """
  function to make a histogram of the number of triggers per time slide
  
  @param slide_trigs: dictionary of time slide triggers
  @param zerolag_trigs: coincInspiralTable
  @param ifos: list of ifos
  """

  nevents = []
  slides = []
  for slide in slide_trigs:
    if ifos:
      nevents.append( slide["triggers"].coinctype(ifos[0],ifos[1]).nevents() )
    else:  
      nevents.append(slide["triggers"].nevents())
    slides.append(slide["slide_num"])
 
    
  hist(nevents)
  figtext(0.13,0.8, " mean = %6.3e" % mean(nevents))
  figtext(0.13,0.75,"sigma = %6.3e" % std(nevents))
  if zerolag_trigs:
    hold(True)
    if ifos:
      nfgevents = zerolag_trigs.coinctype(ifos[0],ifos[1]).nevents()
    else:
      nfgevents = zerolag_trigs.nevents()
    figtext(0.13,0.70,"zero lag = %6.3e" % nfgevents )
    axvline(nfgevents,color='r',linewidth=2)
  
  xlabel('Number of triggers',size='x-large')
  title_text = 'Histogram of number coincident '
  if ifos:
    for ifo in ifos:
      title_text += ifo + ' '
  title_text += 'triggers per time slide'
  title(title_text, size='x-large')
  
  
######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def plotslides(slide_trigs, zerolag_trigs = None, ifos = None):
  """
  function to make a histogram of the number of triggers per time slide
  
  @param slide_trigs: dictionary of time slide triggers
  @param zerolag_trigs: coincInspiralTable
  @param ifos: list of ifos
  """
  nevents = []
  slides = []
  for slide in slide_trigs:
    if ifos:
      nevents.append( slide["triggers"].coinctype(ifos[0],ifos[1]).nevents() )
    else:  
      nevents.append(slide["triggers"].nevents())
    slides.append(slide["slide_num"])
 
  mean_events = mean(nevents)
  std_events = std(nevents)
  plot(slides,nevents,'bx',markersize=12)
  axhline(mean_events,color='k',linewidth=2)
  axhline(mean_events + std_events,color='k',linestyle='--',linewidth=2)
  axhline(mean_events - std_events,color='k',linestyle='--',linewidth=2)
 
  if zerolag_trigs:
    hold(True)
    if ifos:
      nfgevents = zerolag_trigs.coinctype(ifos[0],ifos[1]).nevents()
    else:
      nfgevents = zerolag_trigs.nevents()
    plot([0],[nfgevents],'rx',markersize=12)
 
  
  xlabel('Number of time slide',size='x-large')
  ylabel('Number of triggers in slide',size='x-large')
  title_text = 'Plot of number coincident '
  if ifos:
    for ifo in ifos:
      title_text += ifo + ' '
  title_text += 'triggers per time slide'
  title(title_text, size='x-large')

   

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
