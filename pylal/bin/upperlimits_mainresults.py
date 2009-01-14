#!/usr/bin/python


import sys,os,math
from optparse import *
import ConfigParser
import random
import shutil

##############################################################################
usage = """
usage: %prog [options] 

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )
  parser.add_option("-f", "--config-file",action="store",type="string",\
      metavar=" FILE", help="use configuration file FILE")

  (options,args) = parser.parse_args()
  return options, sys.argv[1:]

def make_dir(directory):
  if not os.path.isdir(directory):
    os.makedirs(directory)

# ============================================================================
# -- get command line arguments
opts, args = parse_command_line()
if not opts.config_file:
  print >> sys.stderr , 'You must specify a config file'
  sys.exit(1)

###################################

cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)
runs = []
spins = []
sourceTypes = []
userTag = cp.get('main','user-tag')
combinedResults = False
if cp.has_option('run-options','run-spin'):
  spins.append('spin')
if cp.has_option('run-options','run-non-spin'):
  spins.append('nonspin')
if cp.has_option('run-options','run-gaussian'):
  for item in cp.options('gaussian-types'):
    sourceTypes.append(item)
if cp.has_option('run-options','run-total-mass'):
  for item in cp.options('total-mass-ranges'):
    sourceTypes.append(item)
if cp.has_option('run-options','run-component-mass'):
  for item in cp.options('component-mass-ranges'):
    sourceTypes.append(item)

for spin in spins:
  for sourceType in sourceTypes:
    dir = spin + '/' + sourceType + '/'
    runs.append((dir,spin,sourceType))

for dir,type,item in runs:
  outdir = 'results/' + dir
  if combinedResults:
    resultDir = 'combine_posteriors/' + dir
  else:
    resultDir = 'plotnumgalaxies_files/' + dir
  resultFile = resultDir + '/' + userTag + '_' + item + '_' + type
  if not combinedResults:
    resultFile += '_combos'
  resultFile += '-combined-posterior.txt'
  outFile = outdir + '/' + userTag + '_posterior-pdf.txt'
  make_dir(outdir)
  shutil.copy(resultFile,outFile)

make_dir('key_plots')
combinedPlotulvsmassDir = 'combine_posteriors/'
combinedPlotulvsmassFiles = ['mcomp_nonspin-combined-rate-v-mass.png']
combinedPlotulvsmassFiles.append('mcomp_spin-combined-rate-v-mass.png')
combinedPlotulvsmassFiles.append('mtotal_nonspin-combined-rate-v-mass.png')
combinedPlotulvsmassFiles.append('mtotal_spin-combined-rate-v-mass.png')
for file in combinedPlotulvsmassFiles:
  shutil.copy(combinedPlotulvsmassDir + '/' + file, 'key_plots/' + file)

combinePosteriorDir = 'combine_posteriors/'
for item in ('bbh','bns','nsbh'):
  for type in ('spin','nonspin'):
    dir = combinePosteriorDir + '/' + type + '/' + item + '/'
    file = userTag + '_' + item + '_' + type + '-posterior-comparison.png'
    shutil.copy(dir +  file, 'key_plots/' + file)

