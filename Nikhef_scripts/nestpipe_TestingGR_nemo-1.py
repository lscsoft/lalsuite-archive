################################################################################

#

#   USER DEFINED VARIABLES

#

################################################################################



xmlname = "injections_dphi6_7000.xml"

parsername = "parser.ini"

basefolder = "/home/salvatore.vitale/IMR_runs/20Hz/IMRPhenom_dphi6/7000"

postprocfolder = "/home/salvatore.vitale/WWW/Testing_GR/20Hz/IMRPhenom_dphi6/7000"

scratchdir = "/scratch/salvatore.vitale/Testing_GR/20Hz/IMRPhenom_dphi6/7000"

logdir = "/people/salvatore.vitale/Testing_GR/20Hz/IMRPhenom_dphi6/7000"

inspnest_dataseed=7000000 ## Don't forget to change that too, increasing by 1000!



hypotheses =["dphi1","dphi2","dphi3","dphi4"]

allPNparams = ["dphi0","dphi1", "dphi2", "dphi3", "dphi4", "dphi5", "dphi6", "dphi7", "dphi8", "dphi9", "spin1z", "spin2z"]



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

events=[0:199]

seed=1

data_seed="""+str(inspnest_dataseed)+"""

analysis-chunk-length=20.0



[condor]

inspnest=/home/salvatore.vitale/opt/lalapps/bin/lalapps_inspnest

combinez=/home/salvatore.vitale/opt/lalapps/bin/lalapps_combine_evidence

datafind=/home/salvatore.vitale/opt/glue/bin/ligo_data_find

mergescript=/home/salvatore.vitale/opt/lalapps/bin/lalapps_merge_nested_sampling_runs

resultspage=/home/salvatore.vitale/opt/pylal/bin/cbcBayesPostProc.py

segfind=/home/salvatore.vitale/opt/glue/bin/ligolw_segment_query

ligolw_print=/home/salvatore.vitale/opt/glue/bin/ligolw_print

coherencetest=/home/salvatore.vitale/opt/lalapps/bin/lalapps_coherence_test



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

	parser_fileobject.write(parser_text)

	parser_fileobject.close()



	os.system("cp "+curdir+"/"+xmlname+" "+basefolder+foldername)



	os.system("lalapps_nest_pipe_reduced -p " + logdir+foldername + " -l " + scratchdir+foldername + " -i "+ basefolder+foldername+"/"+parsername+ " -I "+basefolder+foldername+"/"+xmlname+ " -r " + basefolder+foldername + " --condor-submit")

	#print "lalapps_nest_pipe -i "+ basefolder+foldername+"/"+parsername+ " -I "+basefolder+foldername+"/"+xmlname+ " -r " + basefolder+foldername + " --condor-submit"







# RETURN TO CURRENT WORKING DIRECTORY

os.chdir(curdir)

