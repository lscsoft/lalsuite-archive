
#!/usr/bin/python
import sys
import exceptions
from optparse import *
import glob
import re

from glue import lal
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal import SnglInspiralUtils
from pylal import SimInspiralUtils
from pylal.tools import XLALCalculateEThincaParameter

from pylab import*
import pylal.itertools
import numpy



usage= """
usage: %prog [options]

Takes a file in cache format, sieves it using the pattern if the latter is specified, then writes two new files in cache format: one containing the entries pointing to the files that are physically present and another containing the entries pointing to the files that were not found"""
###############################################################################
# Options to read in Input
###############################################################################
parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )

parser.add_option("-a","--cache-file",action="store",type="string",\
    default=None, metavar="CACHEFILE",help="name of the cache file including the path" )

parser.add_option("-b","--found-cache",action="store",type="string",\
    default=None, metavar="FOUNDFILE",help="name of the cache file in which the found cache entries will be written to" )

parser.add_option("-c","--missed-cache",action="store",type="string",\
    default=None, metavar="MISSEDFILE",help="name of the cache in which the cache entries that were not found will be written to ")

parser.add_option("-d","--sieve-pattern",action="store",type="string",\
    default=None, metavar="SIEVEPATTERN",help="the pattern that will be used to sieve given cache file prior running the check on its entries. If this option is not given, no seiving is performed" )

(opts,args) = parser.parse_args()


if not opts.cache_file:
   print >>sys.stderr, "The option --cache-file must be specified "
   sys.exit(1)

if not opts.found_cache:
   print >>sys.stderr, "The option --found-cache  must be specified "
   sys.exit(1)

if not opts.missed_cache:
   print >>sys.stderr, "The option --missed-cache must be specified "
   sys.exit(1)






AllCache = lal.Cache().read(opts.cache_file)
if opts.sieve_pattern:
  PatternCache = AllCache.sieve(description = opts.sieve_pattern)
else:
  PatternCache = AllCache
PatternCache_Found, PatternCache_Missed = PatternCache.check()
Found_list = []
Found_list = PatternCache_Found.pfnlist()

Missed_list = []
Missed_list = PatternCache_Missed.pfnlist()

PatternCache_Found.write(opts.found_cache)
PatternCache_Missed.write(opts.missed_cache)

