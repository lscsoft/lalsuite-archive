#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesGraceDBinfo
#
#       Copyright 2015
#       Salvatore Vitale <salvatore.vitale@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#===============================================================================
# Preamble
#===============================================================================

#standard library imports
import sys
import os
import numpy as np

from pylal import bayespputils as bppu
from pylal import git_version

__author__="Salvatore Vitale <salvatore.vitale@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

USAGE='''%prog --gid graceDBid --samples posterior_samples.dat

Upload to the pe section of the graceDB page gid information about a lalinference postprocessing run.
For the moment this is maximum a posteriori and stdev for the parameters in the string pars.

'''

def cbcBayesGraceDBinfo(gid=None,samples=None,analysis='LALInference', bcifile=None,bsnfile=None):

  if gid is None or samples is None:
    print "Must provide both a graceDB id and a posterior samples file\n"
    sys.exit(1)

  import ligo.gracedb.rest 
  g=ligo.gracedb.rest.GraceDb()
  samples=os.path.realpath(samples)
  peparser=bppu.PEOutputParser('common')
  commonResultsObj=peparser.parse(open(samples,'r'))
  try:
    pos = bppu.BurstPosterior(commonResultsObj)
    pars=['frequency','quality','hrss']
  except:
    pos = bppu.Posterior(commonResultsObj)
    pars=['mchirp','q','distance']
  strs=[]
  for i in pars:
    if i in pos.names:
      _,which=pos._posMap()
      str="%s parameter estimation: "%analysis
      if i=='hrss':
        str+='%s -- maxP %.3e stdev %.3e'%(i,pos[i].samples[which][0],pos[i].stdev)
      else:
        str+='%s -- maxP %.3f stdev %.3f'%(i,pos[i].samples[which][0],pos[i].stdev)
      strs.append(str)
  if os.path.isfile(bcifile):
    bci=np.loadtxt(bcifile)
  else: bci=None
  if bci is not None:
    str="%s parameter estimation: "%analysis
    str+='logBayes Coherent/Incoherent %.3f'%(bci)
    strs.append(str)

  if os.path.isfile(bsnfile):
    bsn=np.loadtxt(bsnfile)
  else: bsn=None
  if bsn is not None:
    str="%s parameter estimation: "%analysis
    str+='logBayes Signal/Noise %.3f'%(bsn[0])
    strs.append(str)


  for str in strs:
    g.writeLog(gid,str,filename=None,tagname='pe')


if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser(USAGE)
    parser.add_option("-g","--gid", dest="gid",help="GraceDB id", metavar="G123456",default=None)
    parser.add_option("-s","--samples",dest="samples",help="posterior_samples.dat",default=None)
    parser.add_option("--analysis",help="Prefix to use for the graceDB entries. Should be the name of the analysis (default LALInference)",default='LALInference')
    parser.add_option("--bci",dest="bci",help="coherence test file: bci.dat",default=None)
    parser.add_option("--bsn",dest="bsn",help="evidence file: bsn.dat",default=None)
    (opts,args)=parser.parse_args()
    if opts.gid is None:
      print "Must provide a graceDB id with --gid/-g "
      sys.exit(1)
    if opts.samples is None:
      print "Must provide lalinference posterior samples with --samples/-s "
      sys.exit(1)
    cbcBayesGraceDBinfo(opts.gid, opts.samples,analysis=opts.analysis,bcifile=opts.bci,bsnfile=opts.bsn)
