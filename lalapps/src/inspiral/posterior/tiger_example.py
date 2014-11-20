from optparse import OptionParser
import commands
import os
import getpass
import ConfigParser
user=getpass.getuser()

## 2013 Salvatore Vitale. Python to setup both injection and init files for TIGER runs. 
## Will create the injection, the configs, the folders, and call the pipeline exec file.
## 2014 User can pass pre-existing xml table


################################################################################
#
#   USER DEFINED VARIABLES
#
################################################################################

seed=5000          ## Your inspinj seed. The inspnest dataseed will be created from this, adding three zeros at the end (e.g. inspinj 7001 --> inspnest 7001000)
runname="mg_15dot5" # A descriptive name of what is being run. Will be the name of the postprocessing result folder and a part of the injection file, for the records
injfile=None       ## If the user has already an xml table the full path to the table can be passed here, as a string. Otherwise a table will be generated with the pars below

#*** If the user doesn't have an xml table and wants to create one, the following options are used to control what is being injected. If injfile is not None, they are ignored ***#
type_inj="MG"      ## This has to be either GR or MG for modified gravity
shift=15.5         ## This is in percent. If type_inj is GR this will be ignored (you don't need to set it to zero or empty string)
num_events=20 ## This is the number of signals created in the xml file. Inspnest will analize all of them.
# *** #

dic_inj={
#"min_snr":8,
#"max_snr":25, # NetSNR
#"start_freq":30, # f_low for SNR calculation
#"coinc_snr":5.5,
"min_m":1,
"max_m":3,
"seed":seed,
"waveform":'PPEthreePointFivePN',
"min-distance":100000,
"max-distance":150000,
"d-distr":'volume'
}

dic_nest={
"path_to_bin":"/home/%s/lalsuite/tiger_lalinference_tagged/opt/bin/"%user,
}

########### We should not need to change anything below this line, once technical params have been decided

uname=os.uname()
if  any(["atlas" in i for i in uname]):
    WWWdir="WWW/LSC"
else:
    WWWdir="public_html"

if type_inj=="MG":
    dic_inj.update({   
    'enable-mg':'',
    'loglambdaG':shift 
    })

dic_nest.update({
"wwwdir":WWWdir,
"nlive":1024,
"maxmcmc":2048,
"distance-max":150,
"seglen":64,
"srate":4096,
"flow":20,
"user":user,
"dataseed":int(str(dic_inj['seed'])+"000"),  
"runname":runname,
})

logd=None
scrd=None

if  any(["uwm" in i for i in uname]):
    logd=logdir
    scrd=scratchdirlogd=None

   
basefolder=os.getcwd()

#scratchdir = "/scratch/vdbroeck/TaylorT4_TaylorF2/TaylorF2inj/%s/%s"%(type_name,inspinj_seed)            ##
#logdir = "/people/vdbroeck/TaylorT4_TaylorF2/TaylorF2inj/%s/%s"%(type_name,inspinj_seed)                  ## logdir and scratchdir are ignored in all the clusters but UWM.                 

hypotheses =["dchi1","dchi2","dchi3"]

################################################################################
#
#   CREATING THE INJECTION FILE
#
################################################################################

if not injfile:

  print "Creating the xml file\n"
  print dic_nest['runname']
  if not os.path.isdir(basefolder):
      os.makedirs(basefolder)
      
  #inspinjname=os.path.join(basefolder,'injections_%s_%s_SNR_%s_%s.xml'%(dic_nest['runname'],dic_inj['seed'],dic_inj['min_snr'],dic_inj['max_snr']))
  inspinjname=os.path.join(basefolder,'injections_%s_%s.xml'%(dic_nest['runname'],dic_inj['seed']))

  end_time=940000000
  sta_time=930000000

  time_step=((end_time-sta_time)/(num_events-1))

  dic_inj.update({
  "gps-start-time":sta_time,
  "gps-end-time":end_time,
  "time-step":time_step,
  "output": inspinjname,
  "min_mtot":(dic_inj["min_m"]*2),
  "max_mtot":(dic_inj["max_m"]*2)
  })

  #string="lalapps_inspinj --f-lower 10.0 --gps-start-time GPS-START-TIME --gps-end-time GPS-END-TIME --seed SEED --waveform WAVEFORM --l-distr random --i-distr uniform --min-mass2 MIN_M --max-mass2 MAX_M --min-mass1 MIN_M --max-mass1 MAX_M --m-distr componentMass --min-mtotal MIN_MTOT --max-mtotal MAX_MTOT --disable-spin --amp-order 0 --time-step TIME-STEP --disable-spin --output OUTPUT --min-snr MIN_SNR --max-snr MAX_SNR --snr-distr volume  --ligo-fake-psd LALSimAdLIGO --virgo-fake-psd LALSimAdVirgo --ligo-start-freq START_FREQ --virgo-start-freq START_FREQ --ifos H1,L1,V1"
  string="lalapps_inspinj --f-lower 10.0 --gps-start-time GPS-START-TIME --gps-end-time GPS-END-TIME --seed SEED --waveform WAVEFORM --l-distr random --i-distr uniform --min-mass2 MIN_M --max-mass2 MAX_M --min-mass1 MIN_M --max-mass1 MAX_M --m-distr componentMass --min-mtotal MIN_MTOT --max-mtotal MAX_MTOT --disable-spin --amp-order 0 --time-step TIME-STEP --disable-spin --output OUTPUT "

  for p in dic_inj.keys():
      if not p.upper() in string:
          if type(dic_inj[p])==str:
              string=string+ " --"+p +" "+dic_inj[p]
          else:
              string=string+ " --"+p +" "+repr(dic_inj[p])
      else:
          string=string.replace(p.upper()+" ", "%s "%(dic_inj[p]))

  os.system(string)

else:

  dic_inj.update({"output":injfile})

################################################################################
#
#   DEFINE USEFUL FUNCTIONS
#
################################################################################

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def createCombinationsFromList(list):
        outputlist = []
        outputlist.append("GR")
        for L in xrange(len(list)):
                for subset in combinations(list, L+1):
                        """
                        temp = ""
                        for subsetitem in subset:
                                temp += subsetitem
                        """
                        outputlist.append(subset)

        return outputlist

def ensure_dir(f):
    #d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)

def write_pipe_init(dic,stri):
    import string

    for k in dic.keys():
        dk=dic[k]
        if not type(k)==str:
            k=str(k)
        if not type(dk)==str:
            dk=str(dk)
        stri=string.replace(stri,string.upper(k),dk)

    #print stri
    ofile=open("pipeline.ini","w")
    ofile.write(stri)
    ofile.close()

def get_starting_string():
    
    stri="""
[analysis]
ifos=['H1','L1','V1']
engine=lalinferencenest
nparallel=5
dataseed=DATASEED
coherence-test=False

[paths]
#webdir is the base output dir for results pages
webdir=/projects/USER/WWWDIR/RUNNAME
baseurl=https://atlas1.atlas.aei.uni-hannover.de/~USER/WWWDIR/RUNNAME

[input]
# Maximum length to use for automatically-determined psdlength options
max-psd-length=1024
# spacing between trigger times and start of psd estimation
padding=16

timeslides=false
# Uncomment the following line to ignore science segments. Useful when using fake noise
ignore-science-segments=True
events = all

# Section used by the datafind jobs
[datafind]
url-type=file
types={'H1':'H1_LDAS_C02_L2','L1':'L1_LDAS_C02_L2','V1':'HrecOnline'}

[data]
# S5 has LSC-STRAIN, S6 has LDAS-STRAIN
channels={'H1':'H1:LDAS-STRAIN','L1':'L1:LDAS-STRAIN','V1':'V1:h_16384Hz'}

[condor]
lalinferencenest=PATH_TO_BIN/lalinference_nest
lalinferencemcmc=PATH_TO_BIN/lalinference_mcmc
segfind=PATH_TO_BIN/ligolw_segment_query
datafind=PATH_TO_BIN/ligo_data_find
resultspage=PATH_TO_BIN/cbcBayesPostProc.py
ligolw_print=PATH_TO_BIN/ligolw_print
mergescript=PATH_TO_BIN/lalapps_nest2pos
coherencetest=PATH_TO_BIN/lalapps_coherence_test
mpirun=/bin/true
gracedb=/bin/true

[lalinference]
fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}
flow = {'H1':40,'L1':40,'V1':40}

[engine]

# All options in this section are passed to lalinference_nest, lalinference_mcmc,
# and lalinference_bambi. Some useful ones are noted below.
# Options passed to a sampler which does not know them will simply be ignored

# REQUIRED SETTINGS:
# Nlive specifies the number of live points for each job
nlive=NLIVE

# Sampling rate for data
srate=SRATE

# Segment length to use for analysis (should be long enough for whole template
seglen=SEGLEN

# OPTIONAL SETTINGS:

# Use lalinference_nest's resume ability if the run is interrupted
resume=
# lalinference_bambi automatically resumes, use this is you want to force a start from scratch
#noresume=

# approx can be used to manually specify an approximant
# If this is not given, the code will use whatever was injected in the case of a software injection
# Or TaylorF2threePointFivePN if no injection was given.
approx=TaylorF2TestthreePointFivePN

# Control the amplitude order (default: max available)
#amporder=0

# maxmcmc set the maximum chain length for the nested sampling sub-chains. Default 5000
#maxmcmc=5000 # Auto determination is on, but the length cannot be longer than that.

maxmcmc=MAXMCMC

# Number of independent sample to obtain before quitting for lalinference_mcmc
Neff=1000

# Priors
# For all parameters known to lalinference, the min and max default prior can be overwritten with
#parname-min=MIN
#parname-max=MAX

# The starting point for the MCMC chain(s) can be specified with
#parname=START

# Parameters can be fixed to some value with
#fix-parname=FIXVALUE

#currently known parameters, together with default [min-max] are (radiants for angle, Mpc for distance, Msun for masses)

#time                         Waveform time [trigtime-0.1-trigtime+0.1]
#chirpmss                     Chirpmass [1,15.3]
#eta                          Symmetric massratio (needs --use-eta) [0,0.0312]
#q                            Asymmetric massratio (a.k.a. q=m2/m1 with m1>m2) [0.033,1]
#phase                        Coalescence phase [0,2Pi]
#theta_jn                     Angle between J and line of sight [0,Pi]
#distance                     Distance [1,2000]
#logdistance                  Log Distance (requires --use-logdistance) [log(1),log(2000)]
#rightascension               Rightascension [0,2Pi]
#declination                  Declination [-Pi/2,Pi/2]
#polarisation                 Polarisation angle [0,Pi]

#Spin Parameters:
#a1                           Spin1 magnitude (for precessing spins) [0,1]
#a_spin1                      Spin1 magnitude (for aligned spins) [-1,1]
#a2                           Spin2 magnitude (for precessing spins) [0,1]
#a_spin2                      Spin2 magnitude (for aligned spins) [-1,1]
#tilt_spin1                   Angle between spin1 and orbital angular momentum [0,Pi]
#tilt_spin2                   Angle between spin2 and orbital angular momentum  [0, Pi]
#phi_12                       Difference between spins' azimuthal angles [0,2Pi]
#phi_jl                       Difference between total and orbital angular momentum azimuthal angles [0,2Pi]

#Equation of State parameters (requires --use-tidal or --use-tidalT):
#lambda1                      lambda1 [0,3000]
#lambda2                      lambda2 [0,3000]
#lambdaT                      lambdaT [0,3000]
#dLambdaT                     dLambdaT [-500,500]

# Settings for allowed component masses in Solar Masses, with default values
component-max=3.0
component-min=1.0

distance-max=450

# Allowed total masses in Msun (note, used together with component masses, mc,q,eta priors may lead to inconsistencies. Be careful!)
#mtotal-max=35
#mtotal-min=2

# Setting time prior [secs]
#dt=0.1

# The following three options control various marginalised likelihoods. Use at most one.
# Analytically marginalise over phase (only for Newtonian amplitude orders)
#margphi=
# Semi-analytically marginalise over time
#margtime=
# Semi-analytically marginalise over time and phase (only for Newtonian amplitude orders)
#margtimephi=

# By default the code use spins if the choosen approximant supports spin. NOTE that this include TaylorF2, which will be run with aligned spins.
# Several options, here below,  allows the user to choose a particular setting:

#Disable spin for waveforms which would normally have spins
#disable-spin=

# Only allow aligned spins
#aligned-spin=

# Only allow the heavier mass to spin (can be used together with aligned-spin)
#singleSpin=

# Print progress information throughout the run
progress=

# lalinference_bambi allows you to set a target efficiency and tolerance - these are default values
#eff=0.1
#tol=0.5

# Sample in log(distance) for improved sampling
use-logdistance=
# Sample in eta instead than q=m2/m1
#use-eta=

grtest-parameters=GRTESTPARAMETERS

#####################################################################################
[mpi]
# Settings when running with MPI for lalinference_mcmc or lalinference_bambi

# number of CPUs desired and how much memory on each (MB)
machine-count=8
machine-memory=512

#####################################################################################
[resultspage]
# Settings for the results pages (see cbcBayesPostProc.py --help for more)

# Size of the side of the sky bins in degrees
skyres=0.5

# Do no create 2D plots, which take a while to generate
#no2D=

# Send an email linking to the final page
#email=albert.einstein@ligo.org

#####################################################################################
[segfind]
# URL of the segment database
segment-url=https://segdb.ligo.caltech.edu

#####################################################################################
[segments]
# Names of science segments to use
l1-analyze = L1:DMT-SCIENCE:4
h1-analyze = H1:DMT-SCIENCE:4
v1-analyze = V1:ITF_SCIENCEMODE:7
"""
    return stri



################################################################################
#
#   CREATE ALL POSSIBLE COMBINATIONS FROM A LIST
#
################################################################################

curdir = os.getcwd()
foldernames=""
parser_paths=""
allcombinations = createCombinationsFromList(hypotheses)

for run in allcombinations:

    # DETERMINE FOLDER NAME
    foldername = ""
    subhy=""
    if type(run) is not tuple:
        foldername=run
        subhy=run
    else:
        for item in run:
            foldername += item
            subhy+= item +","
        subhy=subhy[:-1]
    ensure_dir(os.path.join(basefolder,foldername))
    os.chdir(os.path.join(basefolder,foldername))
    foldernames+=foldername+' '
    parser_paths+=str(os.path.join(basefolder,foldername,"pipeline.ini"))+" "
    dic_nest.update({
    'runname':os.path.join(dic_nest['runname'],foldername)
    })
    if subhy!='GR':
        dic_nest.update({'GRtestparameters':subhy})
    else:
        dic_nest.update({'GRtestparameters':''})

    st=get_starting_string()
    write_pipe_init(dic_nest,st)
    

#foldernames=foldernames[:-1]
pipestring="%s -I %s -F %s -r %s" %(os.path.join(dic_nest['path_to_bin'],"lalinference_multi_pipe"),dic_inj['output'],foldernames,basefolder)

if logd is not None:
    pipestring+= "-p %s "%logd
if scrd is not None:
    pipestring+= "-l %s "%scrd
pipestring+=" %s "%parser_paths

print pipestring
os.system(pipestring)

# RETURN TO CURRENT WORKING DIRECTORY
os.chdir(curdir)
