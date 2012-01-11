################################################################################
#
#   LOAD RELEVANT MODULES
#
################################################################################

from optparse import OptionParser
import commands
import os

################################################################################
#
#   USER DEFINED VARIABLES
#
################################################################################

inspinj_seed=7000  ## Your inspinj seed. The inspnest dataseed will be created from this, adding three zeros at the end (e.g. inspinj 7001 --> inspnest 7001000)
type_inj="dphi7"   ## This has to be either GR or the name of the test param (e.g. dphi7)
shift=1            ## This is in percent. If type_inj is GR this will be ignored (you don't need to set it to zero or empty string)
number_of_injs=500 ## This is the number of signals created in the xml file. Inspnest will analize all of them.

if type_inj!='GR':
     type_name=type_inj+'_'+repr(shift)+'pc'
else:
     type_name=type_inj

PATH_TO_OPT="/home/salvatore.vitale/lalsuite_branches/nikhef/opt"  ## Path to the opt folder of your installation
parsername = "parser.ini"

basefolder = "/home/salvatore.vitale/IMR_runs/%s/%s"%(type_name,inspinj_seed)                       ##
postprocfolder = "/home/salvatore.vitale/public_html/IMR_runs/%s/%s"%(type_name,inspinj_seed)       ##
scratchdir = "/scratch/salvatore.vitale/IMR_runs/%s/%s"%(type_name,inspinj_seed)                    ##
logdir = "/people/salvatore.vitale/IMR_runs/%s/%s"%(type_name,inspinj_seed)                         ## You only need to change the path leaving the last two levels are they are (i.e. /%s/%s). The code will add the seed and the type of run (GR or tested param +value of the shift)

hypotheses =["dphi1","dphi2","dphi3","dphi4"]
allPNparams = ["dphi0","dphi1", "dphi2", "dphi3", "dphi4", "dphi5", "dphi6", "dphi7", "dphi8", "dphi9", "spin1z", "spin2z"]

################################################################################
#
#   CREATING THE INJECTION FILE
#
################################################################################

print "Creating the xml file\n"

outname='injections_%s_%s.xml'%(type_name,inspinj_seed)

time_step=1.0e+03
gps_start=932170000
gps_end=gps_start+time_step*(number_of_injs ) 

inspinj_command="""lalapps_inspinj \
--output %s \
--f-lower 20.0 \
--gps-start-time %s \
--gps-end-time %s \
--seed %s \
--waveform IMRPhenomFBTestthreePointFivePN \
--min-distance 3.00e+05 \
--max-distance 1.25e+06 \
--d-distr volume \
--l-distr random \
--i-distr uniform \
--min-mass1 5.0 \
--max-mass1 15.0 \
--min-mass2 5.0 \
--max-mass2 15.0 \
--m-distr componentMass \
--min-mtotal 10.0 \
--max-mtotal 30.0 \
--disable-spin \
--time-step %s \
--amp-order 0"""%(outname,gps_start,gps_end,inspinj_seed,time_step)

if type_inj!='GR':
    inspinj_command+=""" \
    --enable-dphi \
    --%s %s"""%(type_inj,repr(shift/100.0)) 

os.system(inspinj_command)
os.system('seed_line=`cat '+outname+' | grep seed`;temp=`echo ${seed_line%\\"*}`;temp2=`echo ${temp##*\\"}`;echo "The file was created with seed $temp2"')
print "\n"
################################################################################
#
#   SETS INSPNEST DATASEED AND XMLNAME
#
################################################################################

inspnest_dataseed=str(inspinj_seed)+"000"
xmlname = outname

################################################################################
#
#   ARGUMENT PARSING
#
################################################################################
"""
parser = OptionParser()
parser.add_option("-i", "--injectionxml", dest="xmlname", help="Injection xml filename", metavar="FILE")
parser.add_option("-b", "--basefolder", dest="basefolder", metavar="FOLDER", help="output folder NS")
parser.add_option("-p", "--postprocfolder", dest="postprocfolder", metavar="FOLDER", help="post processing folder")
parser.add_option("-a", "--parser", dest="parsername", metavar="FILE", help="parser filename (default=parser.ini)", default="parser.ini")

(options, args) = parser.parse_args()

xmlname = options.xmlname
basefolder = options.basefolder
postprocfolder = options.postprocfolder
parsername = options.parsername
"""

if basefolder[-1] != "/":
	basefolder += "/"

if postprocfolder[-1] != "/":
	postprocfolder += "/"

if logdir[-1] != "/":
        logdir += "/"

if scratchdir[-1] != "/":
        scratchdir += "/"

## To do, this can be avoided using os.path.join to join strings forming a path.

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


################################################################################
#
#   CREATE ALL POSSIBLE COMBINATIONS FROM A LIST
#
################################################################################

curdir = os.getcwd()

allcombinations = createCombinationsFromList(hypotheses)

for run in allcombinations:

	# GETTING CORRECT PINPARAMS
	pinparams = ""
	for item in allPNparams:
		if item not in run:
			if pinparams == "":
				pinparams += item
			else:
				pinparams += ","+item

	# DETERMINE FOLDER NAME
	foldername = ""
	for item in run:
		foldername += item

	ensure_dir(basefolder+foldername)
	os.chdir(basefolder+foldername)

	parser_text ="[analysis]"

	parser_text += \
	"""
padding=1
psd-chunk-length=20.0
nlive=1000
nmcmc=100
nparallel=1
ifos=['H1','L1','V1']
events=all
seed=1
data_seed="""+str(inspnest_dataseed)+"""
analysis-chunk-length=20.0

[condor]
inspnest=PATH_TO_OPT/lalapps/bin/lalapps_inspnest
combinez=PATH_TO_OPT/lalapps/bin/lalapps_combine_evidence
datafind=PATH_TO_OPT/glue/bin/ligo_data_find
mergescript=PATH_TO_OPT/lalapps/bin/lalapps_merge_nested_sampling_runs
resultspage=PATH_TO_OPT/pylal/bin/cbcBayesPostProc.py
segfind=PATH_TO_OPT/glue/bin/ligolw_segment_query
ligolw_print=PATH_TO_OPT/glue/bin/ligolw_print
coherencetest=PATH_TO_OPT/lalapps/bin/lalapps_coherence_test

[datafind]
url-type=file

[data]
types=['LALAdLIGO','LALAdLIGO','LALAdVirgo']
channels=['LALAdLIGO','LALAdLIGO','LALAdVirgo']

[inspnest]
dt=0.1
srate=4096.0
approximant=IMRPhenomFBTest

[inspnest_all]
dmin=1.0
dmax=1500.0
deta=0.1
flow=20.0
amporder=0
phaseorder=7
"""
	parser_text +="pinparams=" + pinparams

	parser_text += \
	"""

[results]
skyres=1.0
"""

	parser_text += "basedir="+postprocfolder+foldername
	parser_text += \
	"""

[segfind]
segment-url=https://segdb.ligo.caltech.edu

[segments]
l1-analyze = L1:DMT-SCIENCE:2
h1-analyze = H1:DMT-SCIENCE:2
v1-analyze = V1:ITF_SCIENCEMODE:6

"""

	parser_fileobject = open("parser.ini","w")
        parser_text=parser_text.replace('PATH_TO_OPT',PATH_TO_OPT)
        parser_fileobject.write(parser_text)
	parser_fileobject.close()

	os.system("cp "+curdir+"/"+xmlname+" "+basefolder+foldername)
        os.system("lalapps_nest_pipe_reduced -p " + logdir+foldername + " -l " + scratchdir+foldername + " -i "+ basefolder+foldername+"/"+parsername+ " -I "+basefolder+foldername+"/"+xmlname+ " -r " +
basefolder+foldername + " --condor-submit")
	#print "lalapps_nest_pipe -i "+ basefolder+foldername+"/"+parsername+ " -I "+basefolder+foldername+"/"+xmlname+ " -r " + basefolder+foldername + " --condor-submit"



# RETURN TO CURRENT WORKING DIRECTORY
os.chdir(curdir)
