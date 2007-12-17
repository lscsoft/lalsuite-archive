#!/ldcg/bin/python


# Script for the manual publishing of GEO frames 
#
# 27/04/06 current version simply outputs results to a file rather than 
# publishing to a database.

import os
import re
import sys
import time
import math
import getopt
import popen2
import glob
from types import *
from glue import statedb

# initialise the command line arguments
run = None
version = None
#frsvpath = '/archive/home/chrism/LSCsoft/src/ligotools/packages/FrContrib/bin/FrStateFetcher'
frsvpath = None
#lscdfpath = '/opt/lscsoft/glue/bin/LSCdataFind'
lscdfpath = None
frametype = None
server = None
dbname = None
log = False
startsec = None
startnano = None
nframes = None

# define the regular expression used for frame file parameter extraction from frame filename
filepat = re.compile(r'^[A-Za-z]+\-\w+\-(\d+)\-(\d+)\..+')

# define the command line arguments
shortop = "r:v:r:f:l:t:d:s:q"
longop = [
  "run=",
  "version=",
  "server=",
  "fetcher=",
  "finder=",
  "frametype=",
  "database=",
  "logfile=",
  "start="
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  print >> sys.stderr, \
    "Usage: %s --run [RUN] --version [VER] --logfile [FILE] --frametype [TYPE]" \
      "--server [SER] --fetcher [PROG] --finder [PROG] --database [DBNAME] --start [GPS]" % sys.argv[0]
  sys.exit( 1 )

# parse the command line arguments
for o, a in opts:
  if o in ("-r", "--run"):
    run = a
  elif o in ("-v", "--version"):
    version = int(a)
  elif o in ("-r", "--server"):
    server = a
  elif o in ("-f", "--fetcher"):
    frsvpath = a
  elif o in ("-l", "--finder"):
    lscdfpath = a
  elif o in ("-t", "--frametype"):
    frametype = a
  elif o in ("-d", "--database"):
    dbname = a
  elif o in ("-s", "--logfile"):
    log = a
  elif o in ("-q", "--start"):
    startsec = int(a)
    startnano = int(0)

if not run or not version or not dbname or not log or not server or not frametype:
  print >> sys.stderr, \
    "Usage: %s --run [RUN] --version [VER] --logfile [FILE] --frametype [TYPE] " \
      "--server [SER] --fetcher [PROG] --finder [PROG] --database [DBNAME] --start [GPS]" % sys.argv[0]
  sys.exit(1)

# UNIX epoch is 1970,Jan,1st,00:00:00, GPS epoch is 1980,Jan,6th,00:00:00, difference is 3657 days = 315964800 sec
UNIXtime=time.time()
timenowsec=int(math.floor(UNIXtime-315964800.0))
temptime=time.gmtime(UNIXtime)
datetime="%s/%s/%s %s:%s:%s" %(temptime[2],temptime[1],temptime[0],temptime[3],temptime[4],temptime[5])

# segpat = re.compile(r'^(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(-?\d+)\s([A-Z][0-9]+)\s(.+)')
# open database ?
sdb = statedb.StateSegmentDatabase(run,dbname)

# open logfile for outputting progress
logfile="%s.log" %(log)
if os.path.exists(logfile):
  try: lf=open(logfile, 'a')
  except:
    print >> sys.stderr, "%s *** Error: Unable to open logfile %s for appending, exiting.\n" %(datetime,logfile)
    sys.exit(1)
else:
  try: lf=open(logfile,'w')
  except:
    print >> stderr, "%s *** Error: Unable to create logfile %s for appending, exiting.\n" %(datetime,logfile)
    sys.exit(1)

# open lastlogfile for extracting last published frame times
lastlogfile="%s.log.last" %(log) 
if os.path.exists(lastlogfile):
  try: llf=open(lastlogfile, 'r')
  except:
    lf.write("%s *** Error: Unable to open last logfile %s\n" %(datetime,lastlogfile))
    sys.exit(1)
  # there should only be one line in this file
  for line in llf:
    temp = line.rstrip().split()
  # extract end time of last published file
  startsec=temp[2]
  startnano=temp[3]
# if no last logfile then make one using start time (if start time defined)
elif startsec:
  lf.write("%s *** Last logfile %s does not exist, making last logfile entry %d %d\n" %(datetime,logfile,startsec,startnano))
  try: llf=open(lastlogfile, 'w')
  except:
    lf.write("%s *** Error: Unable to open last logfile %s, exiting.\n" %(datetime,lastlogfile))
    sys.exit(1)
  llf.write("%d %d %d %d %s\n" %(startsec,startnano,startsec,startnano,datetime))
else:
  lf.write("%s *** Error: no last logfile and no start time specified, exiting." %(datetime))
llf.close()


# define temporary data file 
tempdata="%s.tmp" %(log)

# Use LSCdatafind to find data produced since startsec-startnano
lf.write("%s *** Looking for new GEO frames in range (%s-%s)\n" %(datetime,startsec,timenowsec))
try: os.system("%s --server=%s --observatory=G --type=%s --gps-start-time=%s --gps-end-time=%s --url-type=file --match=archive/frames  > %s" %(lscdfpath,server,frametype,startsec,timenowsec,tempdata))
except:
  lf.write("%s *** Error: failed search for GEO frames using %s, exiting.\n" %(datetime,lscdfpath))
  sys.exit(1)

# open temporary file containing all the new frame filenames
lf.write("%s *** Reading new frame file names\n" %(datetime))
if os.path.exists(tempdata):
  try: nff=open(tempdata,'r') 
  except:
    lf.write("%s *** Error: Unable to read from temporary file %s, exiting.\n" %(datetime,tempdata))
    sys.exit(1)
else:
  lf.write("%s *** Error: Temporary file %s does not exist, exiting.\n" %(datetime,tempdata))
  sys.exit(1)

# open results file for the moment (testing storage of results until we can access the segment database)
#resultsfile="%s.results" %(log)
#try: res=open(resultsfile,'w')
#except:
#  lf.write("%s *** Error: Unable to open results file %s, exiting.\n" %(datetime,resultsfile))
#  sys.exit(1)

# loop over each frame file and publish the state vector
for line in nff:
  # extract the frame path from the temporary file
  temp=str(line.rstrip())
  framefile=''
  # we remove the first 16 characters file://localhost
  for i in range(len(temp)-16):
    framefile="%s%s" %(framefile,temp[i+16])

  # extract file information from filename
  sys.stdout.flush()
  lfn = os.path.basename(framefile)
  gpstime = filepat.search(lfn).groups()
  lfn_start = int(gpstime[0])
  lfn_end = lfn_start + int(gpstime[1])

  # register the database ?
  lf.write("%s *** Opening the database %s for times %d - %d\n" %(datetime,dbname,lfn_start,lfn_end))
  try: sdb.register_lfn(lfn,lfn_start,lfn_end)
  except statedb.StateSegmentDatabaseLFNExistsException:
    lf.write("%s *** Warning: LFN %s already exists in database. Trying insertion anyway.\n" %(datetime,lfn))
    #sys.exit(1)

  # Use FrStateFecther to extract the state vector information fronm the file
  lf.write("%s *** Extracting state information from %s\n" %(datetime,framefile))
  sv = popen2.Popen3(frsvpath + " " + framefile)
  sv.tochild.close()
  r = sv.wait()
  if r:
    lf.write("%s *** Error: FrStateFetcher failed: %s, exiting.\n" %(datetime,sv.childerr.readlines()))
    sys.exit(1)
  segdata = sv.fromchild.readlines()
  #print segdata
  del sv

  # extract information from FrStateFetcher output
  for segdatum in segdata:
    segmatch=segdatum.rstrip().split()
    start = int(segmatch[0])
    start_nano = int(segmatch[1])
    end = int(segmatch[2])  
    end_nano = int(segmatch[3])
    sv = int(segmatch[4])
    segnum = 0
    ifo = str(segmatch[7])
 
    # apply bitmask to statevector
    # use only first 6 bits for times after 824828454
    if start > 824828454:
      bitmask = 64
      if math.fmod(sv,bitmask) == 0:
        sv = 0
      else:
        sv = 1
    # use only first 5 bits for times before 824828454    
    else:    
      bitmask = 32
      if math.fmod(sv,bitmask) == 0:
        sv = 0
      else:
        sv = 1
 
    # publish the state vector 
    lf.write("%s *** Publishing state vector to database %s for times %d.%d - %d-%d\n" %(datetime,dbname,start,start_nano,end,end_nano))
    try: sdb.publish_state(ifo,start,start_nano,end,end_nano,version,sv,segnum)
    except statedb.StateSegmentDatabaseSegmentExistsException:
      lf.write("%s *** Error: Unable to publish to the database %s for times %d.%d - %d.%d.\n" %(datetime,dbname,start,start_nano,end,end_nano))
      #sys.exit(1)

    # at present we output to results file
    #res.write("%s %d %d %d %d %d %d %s\n" %(ifo,start,start_nano,end,end_nano,version,sv,segnum))
    last="%d %d %d %d %s\n" %(start,start_nano,end,end_nano,datetime)

  # set flag meaning new frames have been found    
  nframes=1

# close the temporary file
nff.close()

# if new frames published
if nframes:
  # open last published frame logfile
  try: llf=open(lastlogfile,'w')
  except:
    lf.write("%s *** Error: Unable to open %s for writing last published frame times, exiting." %(datetime,lastlogfile))
    sys.exit(1)
  # write to last logfile
  llf.write("%s" %(last))
  llf.close()
else:
  # output null result to log file
  lf.write("%s *** No URLs found\n" %(datetime))
# finally close the logfile
lf.write("%s *** Done, exiting.\n" %(datetime))
lf.close()

# close results file
#res.close()  

# exit
sys.exit(0)

