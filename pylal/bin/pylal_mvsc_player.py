#!/usr/bin/python
"""Calculates and uses MVSC values from cache files.

Currently being written"""

__author__ = 'Tristan Miller <tmiller@caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'
__prog__ = 'pylal_mvsc_player'

##############################################################################
# import modules
import sys,os,random
from optparse import OptionParser
from mvsc_plots import *
from mvsc_htmlwriter import mvsc_html
from matplotlib import pyplot

##############################################################################
# parse options and arguments

parser = OptionParser()
parser.add_option("-r","--run-name", default=False, \
    help="The name of the run (required). Any name is okay, as long as " + \
    "it's the same one used to generate the pat files")
parser.add_option("-t","--gps-time", default=False, \
    help="The gps time of the run (used in printing out the html file)" )
parser.add_option("-S","--stations", default='H1L1', \
    help="Use the given stations (def = H1L1)" )
parser.add_option("-a","--data-sets",help="Choose which data sets for "+\
    "training and validation. Ex: -a BC trains on B, validates on "+\
    "C, and tests on A+D.  Randomized by default." )

parser.add_option("","--zero-lag",default=False,action="store_true", \
                  help="Test on zero-lag triggers in playground.")
parser.add_option("","--open-box",default=False,action="store_true", \
                  help="Test on zero-lag triggers for full data.")
parser.add_option("","--hardware",default=False,action="store_true", \
    help="Test on hardware injection triggers (currently unimplemented).")

parser.add_option("", "--cache-file", default=False, \
    help="Generate new pat files from the given cache file" )
parser.add_option("-p", "--plot-only", action="store_true", default=False, \
    help="Generate plots only, without training or testing again. (Note "+ \
    "that you have to manually reinput -R and -a if you want the html file " \
    "to have the right ones)")
parser.add_option("-d", "--show-plots", action="store_true", default=False, \
    help="Show all plots at the end of program" )
parser.add_option("-O","--enable-output",action="store_true",\
      default=False, \
      help="enable the generation of the html and cache documents")
                  
parser.add_option("-n", default='30', \
    help="number of Bagger training cycles (def=30)" )
parser.add_option("-l", default='4', \
    help="minimal number of entries per tree leaf (def=4)" )
parser.add_option("-c", default='5', \
    help="""criterion for optimization (def=5)                     
                 1 = correctly classified fraction               
		 2 = signal significance s/sqrt(s+b)             
		 3 = purity s/(s+b)                              
		 4 = tagger efficiency Q                         
		 5 = Gini index                               
		 6 = cross-entropy                               
		 7 = 90% Bayesian upper limit with uniform prior 
		 8 = discovery potential 2*(sqrt(s+b)-sqrt(b))   
		 9 = Punzi's sensitivity s/(0.5*nSigma+sqrt(b))""" )
parser.add_option("-g", default='1', \
    help="""per-event loss for (cross-)validation (def=1)          
		 1 - quadratic loss (y-f(x))^2                   
		 2 - exponential loss exp(-y*f(x))               
		 3 - misid fraction""" )
parser.add_option("-s",default=False,help="Number of parameters to sample " +\
                  "for each tree")

parser.add_option("-b", "--balance-data", action="store_true", default=False, \
    help='Train and validate on balanced data which has equal number of ' + \
        'timeslides'+' and injections. This data is much smaller, ' + \
        'because most timeslides are tossed out. (currently unimplemented)' )
parser.add_option("-k", type='int', default=1, \
    help='Train and test multiple times to estimate uncertainty. '+\
         '(under construction)')
parser.add_option("-R", "--seed",action="store",default=False, \
    help="Use a chosen seed. (currently broken to be unreproducible)")

parser.add_option("-o","--output-path",default="html_files", \
      help="Path where the html files and images would be stored. " +\
      "(def = html_files)")
parser.add_option("-q","--data-path",default="Data", \
      help="path where all data files would be stored (def = Data)")
parser.add_option("-w","--pat-path",default="patfiles", \
      help="path where pat files would be stored (def = patfiles)")

opts,args = parser.parse_args()

##############################################################################
#Initialize important values

#LEGEND:
#filecode: name of the data set (ie E14)
#stationcode: name of the stations used (ie H1L1)
#zeropath: filepath of zero-lag data
#escapecode: name of columns which are not to be incorporated in MVSC statistic
#sstr: option s, if it's included
#gstr: option to choose a seed
#letters: list of data sets in order [training,validation,test,none]
#filename: standard unique filename
  #standard is filecode+stationcode+_+letters+_+options
  #(ie E14H1L1_ABC_n30l4c5g1s4)
  #name of file storing forest: filename+.spr
  #name of file with test results: filename+_test.dat
  #name of file with printout: filename+_info
  #name of file with zero-lag results: filename+_playground.dat
  #   or filename+_fulldata.dat
  #name of html file: filename+.html
  #name of images: filename+_+figname+.png
#filepath: path to directory with forests and test results + /filename
#patpath: path to directory with pat files
#valpath,trainpath,testpath: paths to pat files
#htmlfile: mvsc_html object

if opts.run_name:
    filecode = opts.run_name
else:
    print '***Error*** Missing the run-name'
    sys.exit()

stationcode = opts.stations

if opts.open_box:
    #Following code would prompt user about opening box
#    s = raw_input('Are you sure you want to open the box? (yes/no)')
#    if s != 'yes':
#        print 'Aborting...'
#        sys.exit()

    zeropath = 'patfiles/'+ filecode+'/'+stationcode+'setZeroLag_fulldata.pat'
elif opts.zero_lag:
    zeropath = 'patfiles/'+ filecode+'/'+stationcode+ \
               'setZeroLag_playground.pat'
elif opts.hardware:
    zeropath = 'patfiles/' + filecode + '/' + stationcode + \
               'setHardWare.pat'
else:
    zeropath = False

#If k is greater than 1, display only
if opts.k > 1:
    opts.enable_output = False
    opts.show_plots = True

escapecode = "'event_id,L1V1ethinca,H1V1ethinca"

#remove effective snr?
#escapecode += "L1get_effective_snr(),H1get_effective_snr(),"
#escapecode += "H2get_effective_snr(),V1get_effective_snr()"

for i in range(0,len(stationcode),2):
    escapecode += "," + stationcode[i:i+2] + "eff_distance,"
    escapecode += stationcode[i:i+2] + "get_end()"
escapecode += "'"

if opts.s:
    sstr = ' -s ' + opts.s
else:
    sstr = ''

if not opts.seed:
    #random.seed()
    #opts.seed = str(random.randint(0,1048575))
    #print 'Random seed is: ' + opts.seed
    gstr = ''
else:
    gstr = ' -G ' + opts.seed

if opts.data_sets:
    letters = []
    letters.append(opts.data_sets[0])
    letters.append(opts.data_sets[1])
    letters.append('')
    letterset = ['A','B','C','D']
    for l in letterset:
        if (letters[0] != l) & (letters[1] != l):
            letters[2] += l

else:
    if opts.seed:
        random.seed(opts.seed)
    letterset = ['A','B','C','D']
    letters = []
    letters.append(letterset.pop(random.randint(0,3)))
    letters.append(letterset.pop(random.randint(0,2)))
    letters.append(letterset[0]+letterset[1])
    print 'Training on set ' + letters[0] + ', validating on set ' + \
        letters[1] + ', testing on sets ' + letters[2]
    
    opts.data_sets = letters[0]+letters[1]

#standard filenames
filename = filecode + stationcode + '_' + opts.data_sets 
#if opts.balance_data:
#    filename += 'eq'
#else:
#    filename += 'nq'

filename += '_n'+opts.n+'l'+opts.l+'c'+opts.c+'g'+opts.g
if opts.s:
    filename += 's' + opts.s
    
filepath = opts.data_path + '/' + filename
patpath = opts.pat_path

valpath = patpath + '/'+filecode+ '/' + stationcode + 'set' + letters[0] + \
          'Known.pat'
trainpath = patpath + '/'+filecode+'/'+stationcode + 'set' + letters[1] + \
          'Known.pat'
testpath = patpath + '/'+ filecode+'/'+stationcode + 'sets' + letters[2] + \
          '.pat'

#initialize mvsc_html object
if opts.enable_output:
    htmlfile = mvsc_html(opts,filename)

    htmlfile.set_data_sets(letters)

    if not os.path.isdir(opts.output_path):
        os.system('mkdir ' + opts.output_path)
        os.system('mkdir ' + opts.output_path + '/Images')

##############################################################################
#generate pat files from cache

if opts.cache_file:
    time_pat_start = os.times()[4]

    #call pylal_cache_to mvsc.py
    os.system('pylal_cache_to_mvsc.py --cache-file ' + opts.cache_file )

    time_pat_end = os.times()[4]
    print 'Time spent creating pat files: ' + str(time_pat_end-time_pat_start)

    if not os.path.isdir(patpath):
        os.system('mkdir '+patpath)
        
    if os.path.isdir(patpath + '/' + filecode):
        print 'Overwriting previous pat files in ' + filecode
        os.system('rm '+patpath+'/'+filecode+'/*.pat')
    else:
        os.system('mkdir '+patpath+'/' + filecode)
        
    os.system('mv *.pat '+patpath+'/' + filecode)

    print 'Finished generating pat files'

#make a text file with gps time
if opts.gps_time:
    os.system('echo '+opts.gps_time+' > '+patpath+'/'+filecode+'/gpstime')
    
##############################################################################
#call SprBaggerDecisionTreeApp

if not opts.plot_only:
    if not os.path.isdir(opts.data_path):
        os.system('mkdir ' + opts.data_path)
    
    if opts.k > 1:
        for k in range(1,opts.k+1):
            print 'Training run number', k
            filepathk = filepath + 'run' + str(k)
                
            os.system( 'SprBaggerDecisionTreeApp -a 1 -G -i -z '+ \
                escapecode + \
                sstr + ' -n '+opts.n+ \
                ' -l '+opts.l+' -c '+opts.c+' -g '+opts.g+ \
                ' -t ' + valpath+ \
                ' -d 1 -f '+filepathk+'.spr '+ trainpath + \
                ' > '+filepathk+'_info' )
            os.system( 'SprOutputWriterApp -p 10000 -Z '+escapecode+' -a 1 '+ \
                filepathk+ \
                '.spr '+testpath +' '+filepathk+'_test.dat' )
    else:
        time_tree_start = os.times()[4]
        
        if not os.path.isfile(testpath):
            # concatenate files
            temppath1 = patpath + '/'+filecode+'/'+stationcode + 'set' + \
                        letters[2][0] + 'Known.pat'
            temppath2 = patpath + '/'+filecode+'/'+stationcode + 'set' + \
                        letters[2][1] + 'Known.pat'
            os.system("awk 'NF' " + temppath1 + " > " + testpath )
            os.system("awk 'NR > 2' " + temppath2 + " >> " + testpath)
        
        os.system( 'SprBaggerDecisionTreeApp -a 1 -i -z '+escapecode + \
            gstr+ sstr + \
            ' -n '+opts.n+ \
            ' -l '+opts.l+' -c '+opts.c+' -g '+opts.g+ \
            ' -t ' + valpath + \
            ' -d 1 -f '+filepath+'.spr ' + trainpath + \
            ' > '+filepath+'_info' )
        
        os.system( 'SprOutputWriterApp -p 10000 -Z ' + escapecode +' -a 1 '+ \
            filepath+ \
            '.spr ' + testpath +' '+filepath+'_test.dat' )

        print
        print 'SprBaggerDecisionTreeApp -a 1 -i -z '+escapecode + \
            gstr+ sstr + \
            ' -n '+opts.n+ \
            ' -l '+opts.l+' -c '+opts.c+' -g '+opts.g+ \
            ' -t ' + valpath + \
            ' -d 1 -f '+filepath+'.spr ' + trainpath + \
            ' > '+filepath+'_info'
        print
        print  'SprOutputWriterApp -p 10000 -Z ' + escapecode +' -a 1 '+ \
            filepath+ \
            '.spr ' + testpath +' '+filepath+'_test.dat'
        print
        
        if opts.open_box:
            os.system( 'SprOutputWriterApp -p 10000 -Z ' + escapecode + \
                       ' -a 1 '+ filepath+ \
                       '.spr ' + zeropath +' '+filepath+'_fulldata.dat' )
        elif opts.zero_lag:
            os.system( 'SprOutputWriterApp -p 10000 -Z ' + escapecode + \
                       ' -a 1 '+ filepath+ \
                       '.spr ' + zeropath +' '+filepath+'_playground.dat' )
        elif opts.hardware:
            os.system( 'SprOutputWriterApp -p 10000 -Z ' + escapecode + \
                       ' -a 1 '+ filepath+ \
                       '.spr ' + zeropath +' '+filepath+'_hardware.dat' )
        
        time_tree_end = os.times()[4]
        print 'Time spent training and testing trees: ' + \
              str(time_tree_end-time_tree_start)
else:
    if not os.path.isdir(opts.data_path):
        print 'Data files do not exist!  Try again without -p option.'
        sys.exit()

    if not ( os.path.isfile(filepath + '_test.dat') & \
             os.path.isfile(filepath + '_info') & \
             os.path.isfile(filepath + '.spr') ):
        print 'Data files do not exist!  Try again without -p option.'
        sys.exit()

    if (opts.open_box & (not os.path.isfile( filepath+'_fulldata.dat' ) ) ) |\
       (opts.zero_lag & (not os.path.isfile(filepath+'_playground.dat' ))):
        print 'Zero lag files do not exist!  Try again without -p option.'
        sys.exit()

##############################################################################
#Create plots

data = False

if opts.k > 1:
    afar = 1.0/2000
    
    data = []
    fom = []
    for k in range(1,opts.k+1):
        d,cols,cols2 = patread(filepath+'run'+str(k)+'_test.dat',stationcode )
        data.append(d)
        
        f,treesplits,ts_tr,inj_tr,ts_va,inj_va=inforead( \
            filepath + 'run' + str(k) + '_info' )
        fom.append(f)
        
    mid,upper,lower = ROCcompare( data,cols,afar )
    FOMmean(fom)

    print 'True positive rate: ' + str(mid)
    print 'Confidence interval: [' + str(lower) + ',' + str(upper) + ']'
    
elif opts.show_plots | opts.enable_output:
    time_plot_start = os.times()[4]
    
    data,cols,cols2 = patread( filepath + '_test.dat',stationcode )

    ts_trig_ratio = 50
    
    if opts.open_box:
        zerodata,temp1,temp2 = patread(filepath+'_fulldata.dat',stationcode )
    elif opts.zero_lag:
        zerodata,temp1,temp2 = patread(filepath+'_playground.dat',stationcode )
    elif opts.hardware:
        zerodata,temp1,temp2 = patread(filepath+'_hardware.dat',stationcode )
    else:
        zerodata = None

    if zerodata:
        IFANplot(data,cols,zerodata,ts_trig_ratio)
        if opts.enable_output:
            htmlfile.add_figure('IFAN')
    
    mvsc_cutoff = FARplot(data,cols,zerodata,ts_trig_ratio)
    if opts.enable_output:
        htmlfile.add_figure('FAR')
    
    afar,mvsc_cutoff,effsnr_cutoff = \
        fraction_detected(data,cols,mvsc_cutoff=mvsc_cutoff)
    if opts.enable_output:
        htmlfile.add_figure('Frac_vs_SNR')
        htmlfile.set_op_point(afar,mvsc_cutoff,effsnr_cutoff)

    afar,mvsc_cutoff,effsnr_cutoff = \
        fraction_detected(data,cols,mvsc_cutoff=mvsc_cutoff,distance=True)
    if opts.enable_output:
        htmlfile.add_figure('Frac_vs_effdist')
        htmlfile.set_op_point(afar,mvsc_cutoff,effsnr_cutoff)
        
    snr_vs_chisqr(data,cols,afar,zerodata)
    if opts.enable_output:
        htmlfile.add_figure('missed_inj')
        
    lower,upper = ROCplot(data,cols,ts_trig_ratio=ts_trig_ratio)
    if opts.enable_output:
        htmlfile.add_figure('ROC')
        htmlfile.set_efficiency(lower,upper)
    
    mvsc_vs_effsnr(data,cols,mvsc_cutoff,effsnr_cutoff,zerodata)
    if opts.enable_output:
        htmlfile.add_figure('MVSC_vs_effSNR')

    fom,treesplits,ts_tr,inj_tr,ts_va,inj_va = inforead( filepath + '_info' )
    FOMplot(fom)
    if opts.enable_output:
        htmlfile.add_figure('FOM')
        htmlfile.set_treesplits(treesplits)
        htmlfile.set_trigger_info(ts_tr,inj_tr,ts_va,inj_va)
        
        if zerodata:
            htmlfile.set_zeronum(len(zerodata[0]))
                
    #Best 15 events:
    if opts.enable_output:
        if zerodata:
            events = top_events(zerodata,cols,15,stationcode,zeropath)
        else:
            events = top_events(data,cols,15,stationcode,testpath)
    
        htmlfile.set_top_events(events,cols)

    time_plot_end = os.times()[4]
    print 'Time spent creating plots: ' + \
              str(time_plot_end-time_plot_start)

##############################################################################
#Create html page

if opts.enable_output:
    if os.path.isfile(patpath + '/' + filecode + '/gpstime' ):
        try:
            f = open(patpath + '/' + filecode + '/gpstime' )
        except IOError:
            print '***Error!*** Trouble opening file', filename
        else:
            htmlfile.set_gpstime(f.readline())
            f.close()

    htmlfile.write_file()
    
##############################################################################
#Create xml page

pass

##############################################################################
#Show plots

if opts.show_plots:
    pyplot.show()
