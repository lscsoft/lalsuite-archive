#!/usr/bin/python

from optparse import OptionParser
import glob
import sys
import os
parser = OptionParser()
parser.add_option("-g", "--input-glob",action="store",type="string",\
  default=None,metavar="FILES_TO_GLOB",help="glob for these files")
parser.add_option("-o", "--output-file",action="store",type="string",\
  default=None,metavar="OUTPUT FILE",help="output file")
parser.add_option("","--ignore-IFO-times",action="store",type="string",\
  default='', metavar=" USERTAG",\
  help="comma separated list of IFOTYPE_IFOTIME that should not be included in efficiency calculation e.g. H1H2_H1H2,H2L1_H1H2L1. This option will work only with s5 LV search files that are named IFOTIME_IFOTYPE-*xml" )

(options, args) = parser.parse_args()

if not options.input_glob:
  print >> sys.stderr, "You need to specify an input glob."
  sys.exit(1)

if not options.output_file:
  print >> sys.stderr, "You need to specify an output file."
  sys.exit(1)

inputFiles = glob.glob(options.input_glob)
if not inputFiles:
  print >> sys.stderr, "The input-glob returned no files."
  sys.exit(1)
  
# get rid of certain IFO times if necessary
if options.ignore_IFO_times:

  ifo_times_to_ignore = options.ignore_IFO_times.split(",")
  new_inputFiles = []
  for file in inputFiles:
    tmpifotime=file.split("/")[-1].split("_")[0]
    tmpifotype=file.split("/")[-1].split("_")[1].split("-")[0]
    category="_".join([tmpifotype,tmpifotime])
    if not (category in ifo_times_to_ignore):
      new_inputFiles.append(file)

  inputFiles = []
  inputFiles = new_inputFiles


# Once the time from an ifo-time-category has been added,
# the ifo-time will be added to this list.
ifosincluded = []
# The total time analyzed for each ifo-time will be appended to this.
totaltime = 0.0

# Loop over the input *txt files
for file in inputFiles:
  # Split the file name (if it includes a full path) by "/"
  filepieces = file.split("/")
  for item in filepieces:
    # The last item is the actual .txt file name (the full path is now removed)
    if "txt" in item:
      # The first item in the .txt file name is the ifo time.
      ifotime=item.split("_")[0]
  txtfile = open(file)
  for line in txtfile:
    # Look for the "time analysed line..." and pull out the number of seconds
    if "time analysed for triggers" in line:
      words = line.split(" ")
      analyzedtime = words[6].rstrip()
      # If you have not already included this ifo-time-category, then add it.
      if ifotime not in ifosincluded:
        print str(analyzedtime) + " " + str(ifotime)
        # Add the ifo-time-category to the list so that it isn't re-used.
        ifosincluded.append(ifotime)
        totaltime += int(analyzedtime)
  txtfile.close()
print "The total time = " + str(totaltime)

#creat output directories 
output_dir = options.output_file.split("/")[0]
try:os.mkdir(output_dir)
except: pass

for dir in options.output_file.split("/")[1:-1]:
  output_dir = output_dir + "/" + dir
  try: os.mkdir(output_dir)
  except: pass



# Create the output file.
outfile = open(options.output_file,"w")

# Open one of the original input files and use it as a template to write the new output.
origfile = open(inputFiles[0])

for line in origfile:
  # Replace the time analysed value with the new value for the whole analysis.
  if "time analysed for triggers" in line:
    outfile.write("amount of time analysed for triggers " + str(int(totaltime)) + " sec 0 ns\n")
  else:
    outfile.write(line)

origfile.close()
outfile.close()


