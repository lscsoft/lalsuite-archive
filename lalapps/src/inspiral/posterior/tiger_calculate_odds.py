
# Salvatore Vitale <salvatore.vitale@ligo.org> 2015
# Gets a list of hypothesis and BSN files and write the logOdds GR_modGR into file


import ConfigParser
from optparse import OptionParser,OptionValueError
import sys

usage=""" %prog --infiles hyp1_B.txt hyp2_B.txt ... --hyps hyp1 hyp2 ... --out outfile.txt

Reads Bayes factor of the GR and all modGR hypothesis from files, calculates TIGER's ModGR vs GR logOdds ratio (>0 prefers  modGR) and save it into file

"""

import os
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


if __name__=='__main__':

  parser=OptionParser(usage)
  parser.add_option("-o","--out",action="store",type="string",default=None,help="File to save results into",metavar="tiger_odds_0_939939939.0.txt")
  parser.add_option("-i","--infiles",dest="infiles",action="callback", callback=vararg_callback,help="Space separated list of B files, the order must be the same as --hyps list",default=None,metavar="GR_B.txt dchi1_B.txt")
  parser.add_option("-H","--hyps",dest="hyps",action="callback", callback=vararg_callback,help="Space separated list of hypothesis (including GR), the order must be the same as --infiles list",default=None,metavar="GR phi1")
  (opts,args)=parser.parse_args()

  infiles=opts.infiles
  ninfiles=len(infiles)
  hnames=opts.hyps
  nhnames=len(hnames)
  import sys

  if not ninfiles==nhnames:
    print "Must give equal number of hypotheses (--hyps) and B files (--infiles)\n"
    sys.exit(1)
  import numpy as np
  Nhyps=nhnames-1 # Minus one to get rid of the GR hypothesis, which doesn't enter in the calculation of the pre-coefficient
  from numpy import log
  # starting value just depends on Nhyps
  logO_modgr_gr=-log(2.**Nhyps-1.)
  tmp={}

  for (f,h) in zip(infiles,hnames):
    tmp[h]=float(np.loadtxt(f)[0])

  for h in hnames:
    logO_modgr_gr+=(tmp[h]-tmp['GR'])

  try:
    file=open(opts.out,'w')
    file.write("%.3f\n"%logO_modgr_gr)
    file.close()
    print "logOdds ratio modGR vs GR = %f.\n"%logO_modgr_gr
  except:
    print "Could not open file %s to write\n"%opts.out
    sys.exit(1)
