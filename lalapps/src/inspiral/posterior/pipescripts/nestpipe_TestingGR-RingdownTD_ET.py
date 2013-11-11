# This file sets up the tiger pipeline for ringdown signals in ET.
# The script is divided into 3 main sections:
#   1. Setting up injection file, ini files and folder structures (~line 165)
#   2. Setting up the sub files and creating the overall DAG (~line 653)
#   3. Setting up the html page that collects all the individual
#      event html pages into a single table. (~line 1219)
#
# The user only needs to input values in the
# "User defined variables" section (0.1)
#
# Section 3 also creates a cronscript. This script can be run every N minutes
# by adding the line below to the crontab, using the command crontab -e:
#   */N * * * * /pathtotigerpythonscript/cronscripttiger
# Cron will send you an email at every update unless you add
#   MAILTO=""
# to the top of the crontab
# This script will set the permissions of the script to 700



################################################################################
#
#   0.0 LOAD RELEVANT MODULES
#
################################################################################


import ConfigParser
from optparse import OptionParser

from math import ceil

import os
import sys

################################################################################
#
#   0.1 USER DEFINED VARIABLES
#
################################################################################

testoutputon = False
postprocessing = False
createbuggyhtmlpage = False

inspinj_seed=5000   ## Your inspinj seed. The inspnest dataseed will be created from this, adding three zeros at the end (e.g. inspinj 7001 --> inspnest 7001000)

type_inj="dfreq22"       ## This has to be either GR or the name of the test param (e.g. dtau22)

shift=10            ## This is in percent. If type_inj is GR this will be ignored (you don't need to set it to zero or empty string)

number_of_injs=10 ## This is the number of signals created in the xml file. Inspnest will analize all of them.

psd_file="ET_B_data.txt"
psd_file_path="/home/jeroen/tiger_runs"

#The python script to create the overview html page
htmlpythonscript="/home/jeroen/public_html/html_events.py"

ProjectName="TigerRingdownET"
InjectionName="RingdownTDinj"


# Create type_name, used for folder and file names
# typenameextra will be extra info appended like "10Gpc" or "testrun"
typenameextra = "GreatDistance"
if type_inj!='GR':
    type_name=type_inj+'_'+repr(shift)
    type_name+='pc'
else:
    type_name=type_inj

type_name+="_"+typenameextra


## NOTE: Do not change the '%s' format, only what comes before.

PATH_TO_OPT="/home/jeroen/opt_tiger_ringdown"

weburl = "https://ldas-jobs.ligo.caltech.edu/~jmeidam"

basefolder = "/home/jeroen/tiger_runs/%s/%s/%s/%s"%(ProjectName,InjectionName,type_name,inspinj_seed)

publichtmlfolder = "/home/jeroen/public_html"

## logdir and scratchdir are ignored in all the clusters but UWM.
scratchdir = "/scratch/jeroen/%s/%s/%s/%s"%(ProjectName,InjectionName,type_name,inspinj_seed)
logdir     = "/people/jeroen/%s/%s/%s/%s"%(ProjectName,InjectionName,type_name,inspinj_seed)



testingparameters = ["dfreq22","dfreq33","dtau22"]

allmodes = [21,22,33,44]
allparameters = []
for mode in allmodes:
    allparameters.append("tau%d"%mode)
    allparameters.append("freq%d"%mode)


################################################################################
#
#   0.2 Parameters
#
################################################################################

#The mass distribution determines what mass- and massratio values are required
mdistr="finalMass"

#Type of spin distribution
enable_spin = False
spin_distr_flag = "spin-gaussian"

Nlive=512
seglen=5 #s
Nmcmc=200

flow = 10.0 #Hz

#Distances in redshift for the injections:
min_z = 1.5
max_z = 5.0
#Distances in Mpc for the priors (z-range has to be within this range):
min_distance = 10000.0 #Mpc
max_distance = 50000.0 #Mpc

## Component spins of initial binary
min_spin1 = 0.5
max_spin1 = 0.9
min_spin2 = min_spin1
max_spin2 = max_spin1
mean_spin1 = 0.7
stdev_spin1 = 0.2
mean_spin2 = mean_spin1
stdev_spin2 = stdev_spin1

## mass and massratio prior ranges if mdistr="finalMass"
min_rdmass = 500.0 #Msol
max_rdmass = 1000.0 #Msol
q_min = .3
q_max = 1.

## mass prior ranges if mdistr="componentMass"
min_mass1 = 250.0 #Msol
max_mass1 = 500.0 #Msol
min_mass2 = min_mass1
max_mass2 = max_mass1

if(1):
    #calculate min and max component masses from qmin-qmax and Mmin-Mmax.
    max_mass2 = max_rdmass/(q_max+1.0)
    max_mass1 = max_rdmass - max_mass2
    min_mass2 = min_rdmass/(q_max+1.0)
    min_mass1 = min_rdmass - min_mass2

time_step=1000 #s

gps_start=932170000


# ringdown spin prior
min_rdspin = 0.01
max_rdspin = 0.99






################################################################################
#
#   1.0 END OF USER INPUT -- CREATING THE INJECTION FILE
#
################################################################################

postprocfolder = publichtmlfolder+"/%s/%s/%s/%s"%(ProjectName,InjectionName,type_name,inspinj_seed)


#calculate eta from q
eta_min = q_min/((1.0+q_min)*(1.0+q_min))
eta_max = q_max/((1.0+q_max)*(1.0+q_max))

if (testoutputon):
    print "Creating the xml file\n"

if not os.path.isdir(basefolder):

    os.makedirs(basefolder)

outname=os.path.join(basefolder,'injections_%s_%s.xml'%(type_name,inspinj_seed))

gps_end=gps_start+time_step*(number_of_injs )

inspinj_command=(
str(PATH_TO_OPT)+"/bin/lalapps_inspinj \\\n"
"--output "+str(outname)+" \\\n"
"--f-lower 20.0 \\\n"
"--gps-start-time "+str(gps_start)+" \\\n"
"--gps-end-time "+str(gps_end)+" \\\n"
"--seed "+str(inspinj_seed)+" \\\n"
"--waveform RingdownTD \\\n"
"--min-z "+str(min_z)+" \\\n"
"--max-z "+str(max_z)+" \\\n"
"--d-distr covolume \\\n"
"--l-distr random \\\n"
"--i-distr uniform \\\n"
"--amp-order 0 \\\n"
"--time-step "+str(time_step)+" \\\n"
"--m-distr "+str(mdistr)+" \\\n"
)


if (mdistr == "finalMass"):
    inspinj_command+=(
    "--min-rdmass "+str(min_rdmass)+" \\\n"
    "--max-rdmass "+str(max_rdmass)+" \\\n"
    "--min-mratio "+str(q_min)+" \\\n"
    "--max-mratio "+str(q_max)+" \\\n"
    )
    if (enable_spin):
        inspinj_command+=(
        "--enable-spin \\\n"
        "--"+str(spin_distr_flag)+" \\\n"
        "--min-spin1 "+str(min_spin1)+" \\\n"
        "--max-spin1 "+str(max_spin1)+" \\\n"
        "--min-spin2 "+str(min_spin2)+" \\\n"
        "--max-spin2 "+str(max_spin2)+" \\\n"
        "--stdev-spin1 "+str(stdev_spin1)+" \\\n"
        "--mean-spin1 "+str(mean_spin1)+" \\\n"
        "--stdev-spin2 "+str(stdev_spin2)+" \\\n"
        "--mean-spin2 "+str(mean_spin2)
        )
    else:
        inspinj_command+=(
        "--disable-spin"
        )

elif (mdistr == "componentMass"):
    inspinj_command+=(
    "--min-mass1 "+str(min_mass1)+" \\\n"
    "--max-mass1 "+str(max_mass1)+" \\\n"
    "--min-mass2 "+str(min_mass2)+" \\\n"
    "--max-mass2 "+str(max_mass2)+" \\\n"
    )
    if (enable_spin):
        inspinj_command+=(
        "--enable-spin \\\n"
        "--"+str(spin_distr_flag)+" \\\n"
        "--min-spin1 "+str(min_spin1)+" \\\n"
        "--max-spin1 "+str(max_spin1)+" \\\n"
        "--min-spin2 "+str(min_spin2)+" \\\n"
        "--max-spin2 "+str(max_spin2)+" \\\n"
        "--stdev-spin1 "+str(stdev_spin1)+" \\\n"
        "--mean-spin1 "+str(mean_spin1)+" \\\n"
        "--stdev-spin2 "+str(stdev_spin2)+" \\\n"
        "--mean-spin2 "+str(mean_spin2)
        )
    else:
        inspinj_command+=(
        "--disable-spin"
        )
else:
    print "chosen mass distribution not available for ringdown"
    exit(1)





if type_inj!='GR':
    inspinj_command+=""" \\
--%s %.3f"""%( type_inj, shift/100.0 )


if (testoutputon):
    print inspinj_command

os.system(inspinj_command)

os.system('seed_line=`cat '+outname+' | grep seed`;temp=`echo ${seed_line%\\"*}`;temp2=`echo ${temp##*\\"}`;echo "The file was created with seed $temp2"')

if (testoutputon):
    print "\n"

################################################################################
#
#   1.1 SETS LALINFERENCE DATASEED AND XMLNAME
#
################################################################################



inspnest_dataseed=str(inspinj_seed)+"000"

xmlname = outname


if basefolder[-1] != "/":

    basefolder += "/"



if postprocfolder[-1] != "/":

    postprocfolder += "/"





################################################################################
#
#   1.2 DEFINE USEFUL FUNCTIONS
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





################################################################################
#
#   1.3 CREATE ALL POSSIBLE COMBINATIONS FROM A LIST
#
################################################################################


ifos=["E1","E2","E3"]
types_string="{'E1':'interp:%s','E2':'interp:%s','E3':'interp:%s'}"%(psd_file,psd_file,psd_file)
channels_string=types_string
fakecash_string=types_string

ifos_string="['%s','%s','%s']"%(ifos[0],ifos[1],ifos[2])
l1analyze_string = "%s:DMT-SCIENCE:2"%(ifos[0])
h1analyze_string = "%s:DMT-SCIENCE:2"%(ifos[1])
v1analyze_string = "%s:ITF_SCIENCEMODE:6"%(ifos[2])

curdir = os.getcwd()

foldernames=""

parser_paths=""



allcombinations = createCombinationsFromList(testingparameters)



for run in allcombinations:

    # SET CORECT PARAMETERS

    count=0
    if (run != "GR"):
        testparams = ""
        for param in run:
            testparams += "%s"%param
            if (count < len(run)-1):
                testparams += ","
            count+=1
        GRtestparameters = "GRtestparameters=%s"%testparams
    else:
        GRtestparameters = ""


    # DETERMINE FOLDER NAME

    foldername = ""

    for item in run:
        foldername += item

    ensure_dir(basefolder+foldername)

    os.chdir(basefolder+foldername)



    parser_text = (
    "[analysis]\n"
    "ifos="+str(ifos_string)+"\n"
    "engine=lalinferencenest\n"
    "nparallel=1\n"
    "dataseed="+str(inspnest_dataseed)+"\n"
    "coherence-test=False\n"
    "\n"
    "[paths]\n"
    "#webdir is the base output dir for results pages\n"
    "webdir="+str(postprocfolder+foldername)+"\n"
    "basedir="+str(postprocfolder)+"\n"
    "\n"
    "#baseurl is the www address of the above path\n"
    "baseurl=https://ldas-jobs.ligo.caltech.edu/~jmeidam\n"
    "\n"
    "[input]\n"
    "# User-specified length of the psd. if not specified, will be automatically calculated from segment availability\n"
    "# psd-length=256\n"
    "# User-specified psd start time\n"
    "# psd-start-time=\n"
    "\n"
    "# Maximum length to use for automatically-determined psdlength options\n"
    "max-psd-length=1024\n"
    "# spacing between trigger times and start of psd estimation\n"
    "padding=16\n"
    "\n"
    "\n"
    "# Can manually over-ride time limits here\n"
    "#gps-start-time=\n"
    "#gps-end-time=\n"
    "\n"
    "# Can choose a maximum combined false alarm rate when using a pipedown database\n"
    "#max-cfar=0.1\n"
    "\n"
    "# Can manually specify input files here or over-ride on the command line\n"
    "#gps-time-file=\n"
    "#injection-file=\n"
    "#sngl-inspiral-file=\n"
    "#coinc-inspiral-file=\n"
    "#pipedown-db=\n"
    "timeslides=false\n"
    "# Uncomment the following line to ignore science segments. Useful when using fake noise\n"
    "ignore-science-segments=True\n"
    "\n"
    "[datafind]\n"
    "types="+str(types_string)+"\n"
    "\n"
    "[data]\n"
    "# S5 has LSC-STRAIN, S6 has LDAS-STRAIN\n"
    "channels="+str(channels_string)+"\n"
    "\n"
    "# Options for merge script\n"
    "[merge]\n"
    "# Number of posterior samples to generate. If not given will determine automatically\n"
    "#npos=50000\n"
    "\n"
    "[condor]\n"
    "lalinferencenest=PATH_TO_OPT/bin/lalinference_nest\n"
    "lalinferencemcmc=PATH_TO_OPT/bin/lalinference_mcmc\n"
    "segfind=PATH_TO_OPT/bin/ligolw_segment_query\n"
    "datafind=/bin/true\n"
    "resultspage=PATH_TO_OPT/bin/cbcBayesPostProc.py\n"
    "ligolw_print=PATH_TO_OPT/bin/ligolw_print\n"
    "mergescript=PATH_TO_OPT/bin/lalapps_nest2pos\n"
    "coherencetest=PATH_TO_OPT/bin/lalapps_coherence_test\n"
    "mpirun=\n"
    "gracedb=\n"
    "\n"
    "[resultspage]\n"
    "skyres=0.5\n"
    "# Additional options for the results page\n"
    "# --event is set automatically\n"
    "\n"
    "\n"
    "# LALInferenceMCMC options\n"
    "# --lalinfmcmc is set automatically\n"
    "#downsample=1000\n"
    "#deltaLogL=5\n"
    "#fixedBurnin=100000\n"
    "#oldMassConvention\n"
    "\n"
    "# LALInferenceNest options\n"
    "# --Nlive is set automatically from the lalinferencnest section\n"
    "# --ns is set automatically\n"
    "\n"
    "# Send an email when each page is ready (use with caution!)\n"
    "email=jmeidam@nikhef.nl\n"
    "\n"
    "[lalinference]\n"
    "seglen="+str(seglen)+"\n"
    "fake-cache="+str(fakecash_string)+"\n"
    "\n"
    "[lalinferencenest]\n"
    "approx=RingdownTD\n"
    "nlive="+str(Nlive)+"\n"
    "srate=4096\n"
    "progress=\n"
    "mtotalmin="+str(min_rdmass)+"\n"
    "mtotalmax="+str(max_rdmass)+"\n"
    "massmin="+str(min_rdmass)+"\n"
    "massmax="+str(max_rdmass)+"\n"
    "disable-a=\n"
    "a-min="+str(min_rdspin)+"\n"
    "a-max="+str(max_rdspin)+"\n"
    "Dmin="+str(min_distance)+"\n"
    "Dmax="+str(max_distance)+"\n"
    "symMassRatio=\n"
    "eta-min="+str(eta_min)+"\n"
    "eta-max="+str(eta_max)+"\n"
    "template=LALSimRingdown\n"
    "Nmcmc="+str(Nmcmc)+"\n"
    "flow=["+str(flow)+","+str(flow)+","+str(flow)+"]\n"
    "ringdown=\n"
    "modeldomain=time\n"
    ""+str(GRtestparameters)+"\n"
    "\n"
    "[segfind]\n"
    "segment-url=https://segdb.ligo.caltech.edu\n"
    "\n"
    "[segments]\n"
    "l1-analyze = "+str(l1analyze_string)+"\n"
    "h1-analyze = "+str(h1analyze_string)+"\n"
    "v1-analyze = "+str(v1analyze_string)+"\n"
    "\n"
    "[injections]\n"
    "# options to specify software injections\n"
    "#injection-file=/path/to/file.xml\n"
    "\n" )

    parser_fileobject = open("parser_"+foldername+".ini","w")

    parser_text=parser_text.replace('PATH_TO_OPT',PATH_TO_OPT)

    parser_fileobject.write(parser_text)

    parser_fileobject.close()

    foldernames+=str(foldername)+" "

    parser_paths+=str(os.path.join(basefolder,foldername,"parser_"+foldername+".ini"))+" "


    if os.path.exists("%s/%s"%(psd_file_path,psd_file) ):
        os.system("cp %s/%s %s "%(psd_file_path,psd_file,str(os.path.join(basefolder,foldername))))
    else:
        print "psd file note found"
        exit

logd="None"

scrd="None"

#copy the psd data to the basefolder, where the common dag will be
os.system("cp %s/%s %s "%(psd_file_path,psd_file,str(basefolder)))



for i in os.uname():

    if i.find("uwm")!=-1:

        logd=logdir

        scrd=scratchdir

#execute = "python "+path_to_setup_script+"/Setup_TIGER_ETv2.py -i "+ parser_paths + " -I "+outname+ " -r " + basefolder +" -P "+foldernames+" -p " + logd + " -l " + scrd + " -s " + str(gps_start) + " -e " + str(gps_end) + " -N " + str(number_of_injs)

#print "now executing:", execute
#os.system(execute)



###The following variables are required for creating the dags and subfiles
inits=[]
fnames=[]
for ini in parser_paths.split():
    inits.append(ini)
for name in foldernames.split():
    fnames.append(name)

num_of_inits=len(inits)
num_of_fnames=len(fnames)
Nhyp = num_of_fnames

rootdir = basefolder

injfile = outname










################################################################################
################################################################################
################################################################################
################################################################################
#
#   2.0 SETUP DAGS AND SUBS
#
################################################################################
################################################################################
################################################################################
################################################################################

#import uuid
import time

uniquetimestamp = time.strftime("%d%m%Y%H%M%S")


#######################################################################################
## 2.1 Helper functions
#######################################################################################

#copied from nest_multi_parser_reduced.in
def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

#create directory if it does not exist yet
def createdir(path):
  if not ( os.path.exists( path ) ):
    if (testoutputon):
      print "now creating path:"
      print path
    os.mkdir( path )

#######################################################################################
## 2.2 Functions for creating the subfiles
#######################################################################################

#merge_runs.sub
def write_subfile_nest2pos( rundir, logdir, filename, hypname, inifile ):

  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.readfp(open(inifile))
  if cp==None:
    print "cp or opts are none... exiting \n"
    sys.exit(1)

  if (testoutputon):
    print "merge subfile will be created in: "
    print rundir
    print "merge subfile name: "+filename

    if ( cp.has_option('lalinferencenest','GRtestparameters') ):
      print "GRtestparameters check: "+cp.get('lalinferencenest','GRtestparameters')
    else:
      print "GRtestparameters check: None"

  Nlive = cp.get('lalinferencenest','nlive')
  execpath = cp.get('condor','mergescript')

  filetext = (
  "universe = vanilla\n"
  "executable = "+str(execpath)+"\n"
  "arguments = \" --headers $(macroheaders) --Nlive "+str( Nlive )+" --pos $(macropos) $(macroarguments) \"\n"
  "getenv = True\n"
  "log = "+str(logdir)+"/lalinference_pipeline_"+str(hypname)+".log\n"
  "error = "+str(rundir)+"/log/merge-$(cluster)-$(process).err\n"
  "output = "+str(rundir)+"/log/merge-$(cluster)-$(process).out\n"
  "notification = never\n"
  "queue 1\n"
  )

  subfile = open(os.path.join(rundir,filename),"w")
  subfile.write(filetext)
  subfile.close()

#resultspage.sub

def write_subfile_resultspage( rundir, logdir, filename, hypname, inifile ):

  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.readfp(open(inifile))
  if cp==None:
    print "cp or opts are none... exiting \n"
    sys.exit(1)

  if (testoutputon):
    print "resultspage subfile will be created in: "
    print rundir
    print "resultspage subfile name: "+filename

  execpath = cp.get('condor','resultspage')
  skyres = cp.get('resultspage','skyres')
  email = cp.get('resultspage','email')

  filetext = (
  "universe = vanilla\n"
  "executable = "+str(execpath)+"\n"
  "arguments = \" --skyres "+str(skyres)+" --outpath $(macrooutpath) --email "+str(email)+" $(macroarguments) \"\n"
  "RequestMemory = 2000\n"
  "getenv = True\n"
  "log = "+str(logdir)+"/lalinference_pipeline_"+str(hypname)+".log\n"
  "error = "+str(rundir)+"/log/resultspage-$(cluster)-$(process).err\n"
  "output = "+str(rundir)+"/log/resultspage-$(cluster)-$(process).out\n"
  "notification = never\n"
  "queue 1\n"
  )

  subfile = open(os.path.join(rundir,filename),"w")
  subfile.write(filetext)
  subfile.close()

#datafind.sub

def write_subfile_datafind( rundir, logdir, filename, hypname, inifile ):

  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.readfp(open(inifile))
  if cp==None:
    print "cp or opts are none... exiting \n"
    sys.exit(1)

  if (testoutputon):
    print "datafind subfile will be created in: "
    print rundir
    print "datafind subfile name: "+filename

  execpath = cp.get('condor','datafind')

  filetext = (
  "universe = standard\n"
  "executable = "+str(execpath)+"\n"
  "arguments = \" --observatory $(macroobservatory) --url-type file --gps-start-time $(macrogpsstarttime) --gps-end-time $(macrogpsendtime) --output $(macrooutput) --lal-cache --type $(macrotype) \"\n"
  "getenv = True\n"
  "log = "+str(logdir)+"/lalinference_pipeline_"+str(hypname)+".log\n"
  "error = "+str(rundir)+"/log/datafind-$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err\n"
  "output = "+str(rundir)+"/log/datafind-$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out\n"
  "notification = never\n"
  "queue 1\n"
  )

  subfile = open(os.path.join(rundir,filename),"w")
  subfile.write(filetext)
  subfile.close()


#lalinference.sub

def write_subfile_lalinference( rundir, logdir, filename, hypname, inifile):

  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.readfp(open(inifile))
  if cp==None:
    print "cp or opts are none... exiting \n"
    sys.exit(1)

  if (testoutputon):
    print "lalinference subfile will be created in: "
    print rundir
    print "lalinference subfile name: "+filename

  execpath  = cp.get('condor','lalinferencenest')
  approx    = cp.get('lalinferencenest','approx')
  domain    = cp.get('lalinferencenest','modeldomain')
  nlive     = cp.get('lalinferencenest','nlive')
  srate     = cp.get('lalinferencenest','srate')
  eta_min   = cp.get('lalinferencenest','eta-min')
  eta_max   = cp.get('lalinferencenest','eta-max')
  massmin   = cp.get('lalinferencenest','massmin')
  massmax   = cp.get('lalinferencenest','massmax')
  Dmin      = cp.get('lalinferencenest','Dmin')
  Dmax      = cp.get('lalinferencenest','Dmax')
  a_min     = cp.get('lalinferencenest','a-min')
  a_max     = cp.get('lalinferencenest','a-max')
  mtotalmin = cp.get('lalinferencenest','mtotalmin')
  mtotalmax = cp.get('lalinferencenest','mtotalmax')
  nmcmc     = cp.get('lalinferencenest','Nmcmc')
  template  = cp.get('lalinferencenest','template')
  flow      = cp.get('lalinferencenest','flow')
  if ( cp.has_option('lalinferencenest','GRtestparameters') ):
    GRtestparameters_option = "--GRtestparameters "+str(cp.get('lalinferencenest','GRtestparameters') )
  else:
    GRtestparameters_option = ""

  arguments = (
  "--approx "+str(approx)+" "
  "--psdlength $(macropsdlength) "
  "--snrpath "+str(rundir)+"/SNR "
  "--modeldomain "+str(domain)+" "
  "--nlive "+str(nlive)+" "
  "--inj $(macroinj) "
  "--srate "+str(srate)+" "
  "--event $(macroevent) "
  "--seglen $(macroseglen) "
  "--cache $(macrocache) "
  "--trigtime $(macrotrigtime) "
  "--massmin "+str(massmin)+" "
  "--Dmin "+str(Dmin)+" "
  "--template "+str(template)+" "
  "--psdstart $(macropsdstart) "
  "--progress "
  "--Dmax "+str(Dmax)+" "
  "--massmax "+str(massmax)+" "
  "--ringdown "
  "--timeslide $(macrotimeslide) "
  "--mtotalmin "+str(mtotalmin)+" "
  "--Nmcmc "+str(nmcmc)+" "
  "--outfile $(macrooutfile) "
  "--randomseed $(macrorandomseed) "
  "--disable-a "
  "--a-min "+str(a_min)+" "
  "--a-max "+str(a_max)+" "
  "--mtotalmax "+str(mtotalmax)+" "
  "--dataseed $(macrodataseed) "
  "--flow "+str(flow)+" "
  "--symMassRatio "
  "--eta-min "+str(eta_min)+" "
  "--eta-max "+str(eta_max)+" "
  "--channel $(macrochannel) "
  "--ifo $(macroifo) ")

  filetext = (
  "universe = standard\n"
  "executable = "+str(execpath)+"\n"
  "arguments = \" "+str(arguments)+GRtestparameters_option+" \"\n"
  "getenv = True\n"
  "log = "+str(logdir)+"/lalinference_pipeline_"+str(hypname)+".log\n"
  "error = "+str(rundir)+"/log/lalinference-$(cluster)-$(process)-$(node).err\n"
  "output = "+str(rundir)+"/log/lalinference-$(cluster)-$(process)-$(node).out\n"
  "notification = never\n"
  "queue 1\n"
  )

  subfile = open(os.path.join(rundir,filename),"w")
  subfile.write(filetext)
  subfile.close()


def dag_addnode(dagfp, nodename, macros, subfilename, jobpriority = 9999):

  variables = ""

  uniquenodename=uniquetimestamp+nodename

  for macro in macros:
    variables += str(macro[0])+"=\""+str(macro[1])+"\" "

  nodetext = (
  "JOB "+str(uniquenodename)+" "+str(subfilename)+"\n"
  "RETRY "+str(uniquenodename)+" 0\n"
  "VARS "+str(uniquenodename)+" "+str(variables)+"\n"
  )

  if (jobpriority != 9999):

    if(jobpriority > 20 or jobpriority < -20):
      print "error: jobpriority must be <= 20 and >= -20"
      sys.exit(1)

    nodetext += "PRIORITY "+str(uniquenodename)+" "+str(priority)+"\n"

  dagfp.write(nodetext)



#######################################################################################
## 2.3 Create subfiles and dagfiles for each hypothesis
#######################################################################################


for hyp in range(Nhyp):

  inifile = inits[hyp]
  hypname = fnames[hyp]

  if (testoutputon):
    print "========================================"
    print "current hypothesis:",hypname

  sub_nest2pos_name = "merge_runs_"+hypname+".sub"
  sub_datafind_name = "datafind_"+hypname+".sub"
  sub_lalinference_name = "lalinference_"+hypname+".sub"
  sub_resultspage_name = "resultspage_"+hypname+".sub"

  ###
  # Build directory structure for this hypothesis
  ###
  rundir = os.path.join(rootdir, hypname, "rundir")
  logdir = os.path.join(rootdir, hypname, "logdir")
  createdir(rundir)
  createdir(logdir)
  rundir_cashes = os.path.join(rundir, "cashes")
  rundir_engine = os.path.join(rundir, "engine")
  rundir_log = os.path.join(rundir, "log")
  rundir_snr = os.path.join(rundir, "SNR")
  rundir_postsamples = os.path.join(rundir, "posterior_samples")
  createdir(rundir_cashes)
  createdir(rundir_engine)
  createdir(rundir_log)
  createdir(rundir_snr)
  createdir(rundir_postsamples)

  # Also create folder structure for the html script, otherwise it won't properly run the first time:
  createdir(publichtmlfolder)
  createdir(os.path.join(publichtmlfolder,ProjectName))
  createdir(os.path.join(publichtmlfolder,ProjectName,InjectionName))
  createdir(os.path.join(publichtmlfolder,ProjectName,InjectionName,type_name))
  htmlseedpath = os.path.join(publichtmlfolder,ProjectName,InjectionName,type_name,str(inspinj_seed))
  createdir(htmlseedpath)
  createdir(os.path.join(htmlseedpath,hypname))



  write_subfile_nest2pos    ( rundir, logdir, sub_nest2pos_name,     hypname, inifile )
  if(postprocessing): write_subfile_resultspage ( rundir, logdir, sub_resultspage_name,  hypname, inifile )
  write_subfile_datafind    ( rundir, logdir, sub_datafind_name,     hypname, inifile )
  write_subfile_lalinference( rundir, logdir, sub_lalinference_name, hypname, inifile )

  dagname = "lalinference_"+str(hypname)+".dag"
  dagfile = open(os.path.join(rundir, dagname),'w')

  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.readfp(open(inifile))
  #ifos = cp.get('analysis','ifos')
  types = eval(cp.get('datafind','types'))
  ifos = eval(cp.get('analysis','ifos'))


  gpsstart = float(gps_start)
  gpsend   = float(gps_end)
  Nevents  = int(number_of_injs)
  timestep = (gpsend - gpsstart)/Nevents

  # The job priority range in condor, needs to be spread out over all events + a datafind job.
  priority_stepsize = 1
  if (Nevents+1 > 40):
    priority_stepsize = int( ceil(Nevents/40.) )

  priority = 20

  ##########################################
  ## 2.4 Add datafind job to dag for each ifo
  ##########################################



  n=0
  lcffiles = []
  for ifo in ifos:

    lcf = os.path.join(rootdir,"hypname","rundir","cashes",ifo+"-"+str(types[ifo])+"_"+hypname+".lcf")
    lcffiles.append(lcf)

    datafindmacros = [
      ["macrooutput",lcf],
      ["macrogpsendtime",gpsstart],
      ["macroobservatory",ifo],
      ["macrotype",types[ifo]],
      ["macrogpsstarttime",gpsend]
    ]

    nodename = "datafind"+str(n)+"_"+str(hypname)
    dag_addnode(dagfile, nodename, datafindmacros, os.path.join(rootdir,hypname,"rundir",sub_datafind_name), priority )
    #Add a comment for clarity
    dagfile.write("## Job "+nodename+" generates output file "+lcf+"\n")
    n+=1

  dagfile.write("## ------------------------------\n")
  dagfile.write("## End of Datafind jobs\n")
  dagfile.write("## ------------------------------\n")


  for i in range(Nevents):

    if (priority_stepsize > 1):
      if ( i%priority_stepsize == 0 ):
        priority -= 1
    else:
      priority -= 1

    nodename_lalinf = "lalinference_"+str(hypname)+"_event"+str(i)
    nodename_merge  = "mergeruns_"+str(hypname)+"_event"+str(i)
    nodename_results= "resultspage_"+str(hypname)+"_event"+str(i)

    ifo1 = ifos[0]
    ifo2 = ifos[1]
    ifo3 = ifos[2]

    interpstring = "["+str(types[ifo1])+","+str(types[ifo2])+","+str(types[ifo3])+"]"


    #######################
    ## Lalinference node
    #######################

    seed = int(cp.get('analysis','dataseed'))+i
    randomseed = 0 #random seed for the sampling algorithm, will use time if set to 0
    seglen = cp.get('lalinference','seglen')
    trigtime = gpsstart+i*timestep
    psdstart = "932171979"
    psdlength = cp.get('input','max-psd-length')

    lalinf_outfile = os.path.join(rootdir,hypname,"rundir","engine","lalinferencenest_"+str(hypname)+"_event"+str(i)+"-"+ifo1+ifo2+ifo3+"-"+str(trigtime)+".0-"+str(i)+".dat")


    lalinferencemacros = [
    ["macrocache",interpstring],
    ["macroifo","["+str(ifo1)+","+str(ifo2)+","+str(ifo3)+"]"],
    ["macrotimeslide","[0,0,0]"],
    ["macrorandomseed",randomseed],
    ["macrodataseed",seed],
    ["macroseglen",seglen],
    ["macrotrigtime",trigtime],
    ["macrochannel",interpstring],
    ["macropsdstart",psdstart],
    ["macrooutfile",lalinf_outfile],
    ["macroinj",injfile],
    ["macroevent",i],
    ["macropsdlength",psdlength]
    ]

    dag_addnode(dagfile, nodename_lalinf, lalinferencemacros, os.path.join(rootdir,hypname,"rundir",sub_lalinference_name), priority )
    for x in lcffiles:
      dagfile.write("## Job "+nodename_lalinf+" requires input file "+x+"\n")
    dagfile.write("## Job "+nodename_lalinf+" generates output file "+lalinf_outfile+"\n")

    #######################


    #######################
    ## merge_runs node
    #######################

    posterior_filename = os.path.join(rootdir,hypname,"rundir","posterior_samples","posterior_"+ifo1+ifo2+ifo3+"_"+str(trigtime)+"-"+str(i)+".dat")

    mergerunsmacros = [
    ["macropos",posterior_filename],
    ["macroheaders",lalinf_outfile+"_params.txt"],
    ["macroarguments",lalinf_outfile]
    ]

    dag_addnode(dagfile, nodename_merge, mergerunsmacros, os.path.join(rootdir,hypname,"rundir",sub_nest2pos_name), priority )
    dagfile.write("## Job "+nodename_merge+" requires input file "+lalinf_outfile+"\n")
    dagfile.write("## Job "+nodename_merge+" requires input file "+lalinf_outfile+"_params.txt\n")
    dagfile.write("## Job "+nodename_merge+" generates output file "+posterior_filename+"\n")

    #######################


    #######################
    ## resultspage node
    #######################

    webdir = cp.get('paths','webdir')

    resultspage_outpath = os.path.join(webdir,hypname+"_event"+str(i))

    args = (
    posterior_filename+" "
    "--bsn "+posterior_filename+"_B.txt "
    "--inj "+str(injfile)+" "
    "--eventnum "+str(i)
    )

    resultsmacros = [
    ["macrooutpath",resultspage_outpath],
    ["macroarguments",args]
    ]

    if(postprocessing):
        dag_addnode(dagfile, nodename_results, resultsmacros, os.path.join(rootdir,hypname,"rundir",sub_resultspage_name), priority )
        dagfile.write("## Job "+nodename_results+" requires input file "+posterior_filename+"\n")

    #######################


    #######################
    ## CHILD - PARENT
    #######################
    n=0
    for ifo in ifos:
      nodename = "datafind"+str(n)+"_"+str(hypname)
      dagfile.write("PARENT "+uniquetimestamp+nodename+" "+"CHILD "+uniquetimestamp+nodename_lalinf+"\n")
      n+=1
    dagfile.write("PARENT "+uniquetimestamp+nodename_lalinf+" "+"CHILD "+uniquetimestamp+nodename_merge+"\n")
    if(postprocessing): dagfile.write("PARENT "+uniquetimestamp+nodename_merge+" "+"CHILD "+uniquetimestamp+nodename_results+"\n")


    #######################


    dagfile.write("## ------------------------------\n")
    dagfile.write("## End of event"+str(i)+"\n")
    dagfile.write("## ------------------------------\n")



  dagfile.write("## ==============================\n")
  dagfile.write("## End of hypothesis "+str(hypname)+"\n")
  dagfile.write("## ==============================\n")


  dagfile.close()




#####################################################################
### It should now be a simple matter of appending all individual dags
### into a single dag. Thanks to the job priorities, each event is
### handled with the same priority across all hypotheses.
#####################################################################

commondag = os.path.join(rootdir,"common_dag.dag")

if(os.path.isfile(commondag)):
  os.system("rm "+commondag)

for hyp in range(Nhyp):

  hypname = fnames[hyp]
  currentdag = os.path.join(rootdir, hypname, "rundir", "lalinference_"+str(hypname)+".dag")

  if (testoutputon):
    print "adding file "+currentdag+" to common_dag.dag"


  os.system("cat "+currentdag+" >> "+commondag)

print ""
print "TIGER was set up in:"
print basefolder
print "for hypotheses:"
for i in range(Nhyp):
    print fnames[i]
print ""











################################################################################
################################################################################
################################################################################
################################################################################
#
#   3.0 SETUP HTML PAGE
#
################################################################################
################################################################################
################################################################################
################################################################################

from os.path import dirname



if (createbuggyhtmlpage):
    projecthtmlpath = dirname( dirname( os.path.normpath(postprocfolder) ) )
    outputfile = publichtmlfolder+"/%sProject.html"%ProjectName
    projectlink = weburl+"/%s/%s/"%(ProjectName,InjectionName)
    rundir = dirname( dirname( os.path.normpath(basefolder) ) )

    html_executable = "python "+htmlpythonscript+" -p "+projecthtmlpath+" -r "+rundir+" -o "+outputfile+" -u "+projectlink+" -s "+str(gps_start)+" -t "+str(time_step)
    if(postprocessing): html_executable+="-n"

    if (testoutputon):
        print html_executable

    os.system(html_executable)

    with open(os.path.join(dirname(htmlpythonscript),"cronscripttiger"),"w") as f:
        f.write("#!/bin/bash\n#\n")
        f.write(html_executable+"\n")
        f.write("echo \"Updated page on: $(date)\" >> "+os.path.join(dirname(htmlpythonscript),"cronscripttiger.log"))

    chmod700 = "chmod 700 "+os.path.join(dirname(htmlpythonscript),"cronscripttiger")
    os.system(chmod700)

    if (testoutputon):
        print "changed cronscript file permissions:\n"+"chmod 700 "+os.path.join(dirname(htmlpythonscript),"cronscripttiger")

    print "The html project page will be written to:"
    print "%s"%outputfile
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

## Need to run the dag from that folder, otherwise it can't find the psd file
curdir = os.getcwd()
if (os.path.normpath(basefolder) == curdir):
    print "now run condor_submit_dag common_dag.dag"
else:
    print "now run \ncd "+str(basefolder)+" && condor_submit_dag common_dag.dag"
    if (testoutputon): print ""
    if (testoutputon): print "also... Don't forget to source the correct environment!"






