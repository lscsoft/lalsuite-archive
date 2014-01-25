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
"dmax":150,
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
basedir=/projects/USER/WWWDIR/RUNNAME
baseurl=

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
types={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdLIGO'}

[data]
# S5 has LSC-STRAIN, S6 has LDAS-STRAIN
channels={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}

# Options for merge script
[merge]
# Number of posterior samples to generate. If not given will determine automatically
#npos=50000

[condor]
lalinferencenest=PATH_TO_BIN/lalinference_nest
lalinferencemcmc=PATH_TO_BIN/lalinference_mcmc
segfind=PATH_TO_BIN/ligolw_segment_query
datafind=/bin/true
resultspage=PATH_TO_BIN/cbcBayesPostProc.py
ligolw_print=PATH_TO_BIN/ligolw_print
mergescript=PATH_TO_BIN/lalapps_nest2pos
coherencetest=PATH_TO_BIN/lalapps_coherence_test
mpirun=
gracedb=

[resultspage]
skyres=1.0
basedir=/home/USER/WWWDIR/RUNNAME

[lalinference]
seglen=SEGLEN
fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}

[engine]
nlive=NLIVE
srate=SRATE
template=LALSim
progress=
maxmcmc=MAXMCMC
lalsimulationinjection=
flow=[FLOW,FLOW,FLOW]
Dmax=DMAX
dt=0.1
#spinOrder=0
#tidalOrder=0
#inj-tidalOrder=0
#inj-spinOrder=0
amporder=0
approx=TaylorF2TestthreePointFivePN
proposal-no-extrinsicparam=
resume=
margphi=
GRtestparameters=GRTESTPARAMETERS
#pinparams=PINPARAMS
downsample=1000
deltalogl=5

[segfind]
segment-url=https://segdb.ligo.caltech.edu

[segments]
l1-analyze = L1:DMT-SCIENCE:2
h1-analyze = H1:DMT-SCIENCE:2
v1-analyze = V1:ITF_SCIENCEMODE:6
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
