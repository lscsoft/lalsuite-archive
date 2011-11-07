#!/usr/bin/python
__author__ = "Hugo Marrochio <hcmarroc@mit.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

__prog__="auxmvc_plot_mvsc_channels_significance.py"
__Id__ = "$Id$"


import re
import numpy
import readline
import matplotlib
matplotlib.use('Agg')
import pylab
import pdb
import numpy
from pylal import git_version
from pylal import auxmvc_utils
from pylal import InspiralUtils
import bisect
import pickle
import glob
from optparse import OptionParser
import os

# Define functions

# The function calculates the average and rms mean of significance for each channel.
# The max and min values for each channel are stored, as well as the normalized deviation from rms mean
def characterize_channel(channels, significance, n):
    list_channels=[]
    max_variable=[]
    min_variable=[]
    max_dif=[]
    min_dif=[]
    values=[]
    mean_rms=[]
    n_channels=len(channels)/n
    for i in range(n_channels):
        list_channels.append([channels[n*i][:-7],(significance[n*i] +significance[n*i+1] + significance[n*i+2] + significance[n*i+3] + significance[n*i+4])])
        list_channels[i][1]=list_channels[i][1]/(n)
        values.append(list_channels[i][1])
        mean_rms.append(((significance[n*i]**2 +significance[n*i+1]**2 + significance[n*i+2]**2 + significance[n*i+3]**2 + significance[n*i+4]**2)/(n))**0.5)
        max_variable.append(max(significance[n*i],significance[n*i+1],significance[n*i+2],significance[n*i+3],significance[n*i+4]))
        min_variable.append(min(significance[n*i],significance[n*i+1],significance[n*i+2],significance[n*i+3],significance[n*i+4]))
        if list_channels[i][1] ==0:
            max_dif.append(max_variable[i])
            min_dif.append(min_variable[i])
        else:
            max_dif.append((max_variable[i]-mean_rms[i])/mean_rms[i])
            min_dif.append((min_variable[i]-mean_rms[i])/mean_rms[i])
    return list_channels, max_variable, min_variable, max_dif, min_dif, values, mean_rms

# Function that returns a list of which parameter in the channel corresponds to a maximum value (number between 0-4)
def characterize_max(significance, threshold):
    characterize_max=[]
    maxi=[]
    k=0
    print(significance)
    n_channels=len(significance)/5
    for i in range(n_channels):
        print(i, k)
        average=(significance[5*i] + significance[5*i+1] +significance[5*i+2] +significance[5*i+3] +significance[5*i+4])/5
        print (average, threshold)
        if average > threshold:
            characterize_max.append(0)
            maxi.append(significance[5*i])
            dummy=significance[5*i]
            for j in range(4):
                if significance[5*i+j+1] > dummy:
                    characterize_max[k]=j+1
                    maxi[k] = significance[5*i+j+1]
                    dummy = significance[5*i+j+1]
            k=k+1
    return characterize_max

# Function that returns a list of which parameter in the channel corresponds to a minimum value (number between 0-4)
def characterize_min(significance, threshold):
    characterize_min=[]
    mini=[]
    k=0
    print(significance)
    n_channels=len(significance)/5
    for i in range(n_channels):
        print(i, k)
        average=(significance[5*i] + significance[5*i+1] +significance[5*i+2] +significance[5*i+3] +significance[5*i+4])/5
        print (average, threshold)
        if average > threshold:
            characterize_min.append(0)
            mini.append(significance[5*i])
            dummy=significance[5*i]
            for j in range(4):
                if significance[5*i+j+1] < dummy:
                    characterize_min[k]=j+1
                    mini[k] = significance[5*i+j+1]
                    dummy = significance[5*i+j+1]
            k=k+1
    return characterize_min

#Tries to calculate correlation and covariance for two different parameters i and j (between 0-4)
#Not sure of the relevance, just trying to calculate different things.
def cov_matrix(dvals, i, j):
    tempi=[]
    tempj=[]
    p=len(dvals)/5
    for k in range(p):
        tempi.append(dvals[5*k+i])
        tempj.append(dvals[5*k+j])
    X = numpy.vstack((tempi,tempj))
    return numpy.cov(X), numpy.corrcoef(X)
    
    
        
# Print a list of the channel name, as well as its average significance
def print_output(dtest_FOM,dtest,n,thr,name, output_path):
    g= open(output_path + "/"+"list_output_" + name  + ".txt", "w")
    f= open(output_path + "/"+"list_output_" + name  + "_FOM.txt", "w")
    j=-1
    while j >= -(n) and dtest[j][1]>=thr:
        g.write(str(dtest[j][0]) + '  =  ' + str(dtest[j][1]) +'\n')
        f.write(str(dtest_FOM[j][0]) + '  =  ' + str(dtest_FOM[j][1]) +'\n')
        j=j-1
    g.close()
    f.close()
    return
    





parser=OptionParser(usage="Plot histograms of significance rank of auxiliary channels assigned to them by MVSC", version = git_version.verbose_msg)
parser.add_option("-i","--input", default="0" ,help="Provide the path for the input files and globbing pattern.")
parser.add_option("","--histograms", action="store_true", default=False, help="Make Cumulative H
istograms, Histograms and  relative distribution of the channels")
parser.add_option("-n", "--N-most-singificant-channels",type="int", default=None, help="It outputs the n channels with largest significance. As default, all channels will be in the output")
parser.add_option("","--threshold", default ="0.0",type="float", help="The output contains all channels with bigger significance than this specified threshold. Generally used if the user knows the threshold value for significance")
parser.add_option("","--threshold-maximum-minimum", default ="0.0",type="float", help="Threshold used in determining the max and min distribution in the variables")
parser.add_option("","--plotdeviation", action="store_true", default=False, help="Make plots to
test the distance of each variable to the average of the five, as a function of significance")
parser.add_option("-u","--usertag",action="store",type="string", default='',metavar=" USERTAG",help="a user tag for the output file names. Try to specify the name of a particular run")
parser.add_option("-O","--enable-output",action="store_true",default="True",  metavar="OUTPUT",help="enable the generation of the html and cache documents")	 
parser.add_option("", "--figure-resolution",action="store",type="int",default=50, help="dpi of t:he thumbnails (50 by default)")
parser.add_option("", "--html-for-cbcweb",action="store",default=False, metavar = "CVS DIRECTORY", help="publish the html "\
      "output in a format that can be directly published on the cbc webpage "\
      "or in CVS. This only works IF --enable-output is also specified. The "\
      "argument should be the cvs directory where the html file will be placed "\
      "Example: --html-for-cbcweb protected/projects/s5/yourprojectdir")
parser.add_option("","--verbose", action="store_true",\
      default=False, help="print information" )
parser.add_option("-P","--output-path",action="store",\
      type="string",default="",  metavar="PATH",\
      help="path where the plots and will be stored")
#parser.add_option("","--output_path",action="store",\
#      type="string",default="",  metavar="PATH",\
#      help="path where the results will be stored. Generally prefer to set the same as the previous option, for the plots and html file")
(opts,args)=parser.parse_args()

dic={}
dkeys=[]
dval_FOM=[]
dval_splits=[]
dsort=[]
dev_sig=[]
dev_dt=[]
dev_dur=[]
dev_freq=[]
dev_npts=[]
index=0
i=0
list_of_inputs=glob.glob(opts.input)

# First part of the code: Read all files from round-robin and store the significance, average by the number of files read.

for inp in list_of_inputs:
    f=open(inp,'r')
    for line in f:
        words = re.findall(r'\w+', line)
        if words[0] == 'Variable':
            if index == 0:
                k=(words[1]+'-'+words[2])
                v_FOM = float(words[7]+'.'+words[8])/len(list_of_inputs)
                v_splits = float(words[4])/len(list_of_inputs)
                dval_FOM.append(v_FOM)
                dval_splits.append(v_splits)
                dkeys.append(k)
                dsort.append([k,v_FOM])
                i=i+1
            else:
                dval_FOM[i] = dval_FOM[i] + float(words[7]+'.'+words[8])/len(list_of_inputs)
                dval_splits[i] = dval_splits[i] + float(words[4])/len(list_of_inputs)
                dsort[i][1] = dsort[i][1]+float(words[7]+'.'+words[8])/len(list_of_inputs)
                i=i+1
    index=index+1
    i=0
list_channels_FOM, max_variable, min_variable, max_dif, min_dif, values_FOM, mean_rms = characterize_channel(dkeys,dval_FOM,5)
list_channels, max_variable_s, min_variable_s, max_dif_s, min_dif_s, values_splits, mean_rms_s = characterize_channel(dkeys,dval_splits,5)
list_channels_FOM.sort(key=lambda a: a[1])
list_channels.sort(key=lambda a: a[1])
if opts.N_most_significant_channels == None:
    n=len(values_FOM)
else:
    n=opts.N_most_significant_channels
print(n)
print_output(list_channels_FOM, list_channels, n, opts.threshold, opts.usertag, opts.output_path)

max_channel=characterize_max(dval_FOM,opts.threshold_maximum_minimum)
min_channel=characterize_min(dval_FOM, opts.threshold_maximum_minimum)
max_channel_s=characterize_max(dval_splits,opts.threshold_maximum_minimum)
min_channel_s=characterize_min(dval_splits, opts.threshold_maximum_minimum)


X,Y=cov_matrix(dval_FOM, 1, 4)
print(X)
print(Y)

print('next')
X,Y=cov_matrix(dval_FOM, 0, 1)
print(X)
print(Y)




# Initializing the html output
InspiralUtils.message(opts, "Initialisation...")
opts = InspiralUtils.initialise(opts, __prog__, __version__)
fnameList = []
tagList = []
fig_num = 0
comments = ""

# Initialize the plotting options
if opts.plotdeviation == True:

    # 1 
    fig_num += 1
    
    pylab.figure(1)
    pylab.plot(values_FOM, max_dif,'yo')
    pylab.title('Deviation from average for max value using Delta FOM')
    pylab.xlabel('Delta FOM')
    pylab.ylabel('Deviation from average for max value for each channel')
    pylab.savefig(opts.usertag + '_max_dif_FOM.png')

    name = "max_dif_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("max_dif_FOM")

    pylab.close()


    # 2 
    fig_num += 1
    
    pylab.figure(1)
    pylab.plot(values_splits, max_dif_s,'yo')
    pylab.title('Deviation from average for max value using number of splittings')
    pylab.xlabel('Number of splittings')
    pylab.ylabel('Deviation from average for max value for each channel')
    pylab.savefig(opts.usertag + '_max_dif.png')

    name = "max_dif" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("max_dif")

    pylab.close()

    # 3 
    fig_num += 1
    pylab.figure(1)
    pylab.plot(values_FOM,min_dif,'yo')
    pylab.title('Deviation from average for min value using delta FOM')
    pylab.xlabel('Delta FOM')
    pylab.ylabel('Deviation from average for min value for each channel')
    pylab.savefig(opts.usertag + '_min_dif_FOM.png')

    name = "min_dif_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("min_dif_FOM")
    
    pylab.close()

    # 4 
    fig_num += 1
    pylab.figure(1)
    pylab.plot(values_splits,min_dif_s,'yo')
    pylab.title('Deviation from average for min value using number of splittings')
    pylab.xlabel('Number of splittings')
    pylab.ylabel('Deviation from average for min value for each channel')
    pylab.savefig(opts.usertag + '_min_dif.png')

    name = "min_dif" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("min_dif")
    
    pylab.close()



if opts.histograms == True:
    # 5
    fig_num += 1
    pylab.figure(1)
    pylab.hist(values_FOM,150, histtype='stepfilled',cumulative = -1)
    pylab.title('Cumulative histogram for delta FOM')
    pylab.xlabel('Delta FOM')
    pylab.ylabel('Cumulative histogram')
    pylab.savefig(opts.usertag + '_hist_cumulative_FOM.png')

    name = "hist_cumulative_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_cumulative_FOM")
    
    pylab.close()

    # 6
    fig_num += 1
    pylab.figure(1)
    pylab.hist(values_splits,150, histtype='stepfilled',cumulative = -1)
    pylab.title('Cumulative histogram for number of splittings')
    pylab.xlabel('Delta FOM')
    pylab.ylabel('Cumulative histogram')
    pylab.savefig(opts.usertag + '_hist_cumulative.png')

    name = "hist_cumulative" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_cumulative")
    
    pylab.close()


    # 7
    fig_num += 1
    pylab.figure(1)
    pylab.hist(values_FOM,30, histtype='stepfilled')
    pylab.title('Histogram for delta FOM')
    pylab.xlabel('Delta FOM')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist_FOM.png')

    name = "hist_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_FOM")
    
    pylab.close()


    # 8
    fig_num += 1
    pylab.figure(1)
    pylab.hist(values_splits,30, histtype='stepfilled')
    pylab.title('Histogram for number of splittings')
    pylab.xlabel('Number of splittings')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist.png')

    name = "hist" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist")
    
    pylab.close()

    # 9
    fig_num += 1
    pylab.figure(1)
    pylab.hist(max_dif,50, histtype='stepfilled')
    pylab.title('Deviation from the maximum rms average for delta FOM')
    pylab.xlabel('deviation for max-rms_avg for delta FOM')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist_max_FOM.png')

    name = "hist_max_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_max_FOM")
    
    pylab.close()

    # 10
    fig_num += 1
    pylab.figure(1)
    pylab.hist(max_dif_s,50, histtype='stepfilled')
    pylab.title('Deviation from the maximum rms average for number of splittings')
    pylab.xlabel('deviation for max-rms_avg')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist_max.png')

    name = "hist_max" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_max")
    
    pylab.close()
    
    # 11
    fig_num += 1
    pylab.figure(1)
    pylab.hist(min_dif, 50, histtype='stepfilled')
    pylab.title('Deviation from the minimum rms average for delta FOM')
    pylab.xlabel('deviation for min-rms_avg')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist_min_FOM.png')

    name = "hist_min_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_min_FOM")
    pylab.close()

    # 12
    fig_num += 1
    pylab.figure(1)
    pylab.hist(min_dif_s, 50, histtype='stepfilled')
    pylab.title('Deviation from the minimum rms average for number of splittings')
    pylab.xlabel('deviation for min-rms_avg')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_hist_min.png')

    name = "hist_min" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("hist_min")
    pylab.close()


    # 13
    fig_num += 1
    pylab.figure(1)
    pylab.hist(max_channel, bins=(0,1,2,3,4,5), histtype='stepfilled')
    pylab.title('Maximum channel histogram. Each number represents each variable, in the order of the input, for delta FOM')
    pylab.xlabel('Channel for max')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_max_channel_FOM.png')

    name = "max_channel_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("max_channel_FOM")
    
    pylab.close()

    # 14
    fig_num += 1
    pylab.figure(1)
    pylab.hist(max_channel_s, bins=(0,1,2,3,4,5), histtype='stepfilled')
    pylab.title('Maximum channel histogram. Each number represents each variable, in the order of the input, for number of splittings')
    pylab.xlabel('Channel for max')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_max_channel.png')

    name = "max_channel" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("max_channel")
    
    pylab.close()

    # 15
    fig_num += 1
    pylab.figure(1)
    pylab.hist(min_channel, bins = (0,1,2,3,4,5), histtype='stepfilled')
    pylab.title('Minimum channel histogram. Each number represents each variable, in the order of the input, for delta FOM')
    pylab.xlabel('Channel for min')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_min_channel_FOM.png')

    name = "min_channel_FOM" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("min_channel_FOM")
    
    pylab.close()

    # 16
    fig_num += 1
    pylab.figure(1)
    pylab.hist(min_channel_s, bins = (0,1,2,3,4,5), histtype='stepfilled')
    pylab.title('Minimum channel histogram. Each number represents each variable, in the order of the input, for number of splittings')
    pylab.xlabel('Channel for min')
    pylab.ylabel('Histogram')
    pylab.savefig(opts.usertag + '_min_channel.png')

    name = "min_channel" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("min_channel")
    
    pylab.close()


    # 17
    fig_num += 1
    pylab.figure(1)
   
    pylab.plot(values_FOM,values_splits,'yo')
    pylab.title('Correlation between delta FOM and number of splittings')
    pylab.xlabel('Values for delta FOM')
    pylab.ylabel('Number of splits')
    pylab.savefig(opts.usertag + '_FOM_splits.png')

    name = "FOM_splits" 
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append("FOM_splits")
    
    pylab.close()

##############################################################################################################


html_filename = InspiralUtils.write_html_output(opts, args, fnameList, tagList, comment=comments)
InspiralUtils.write_cache_output(opts, html_filename, fnameList)

if opts.html_for_cbcweb:
  html_filename_publish = InspiralUtils.wrifacte_html_output(opts, args, fnameList, tagList, cbcweb=True)

##############################################################################################################
