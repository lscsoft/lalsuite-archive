#!/usr/bin/python


import sys,os,math
from optparse import *
import ConfigParser
import random

##############################################################################
usage = """
usage: %prog [options] 

"""

def determine_reduced_combos(combo):
  cp2 = ConfigParser.ConfigParser()
  cp2.addsection('input')
  ifos = []
  temp = 0
  while temp < len(combo):
    ifos.append(combo[temp:temp+2])
    temp = temp+2
  if ('H1') in ifos:
    cp2.set('input','h1-data','')
  if ('H2') in ifos:
    cp2.set('input','h2-data','')
  if ('L1') in ifos:
    cp2.set('input','l1-data','')
  if ('V1') in ifos:
    cp2.set('input','v1-data','')
  if len(combo) > 5:
    cp2.set('input','two-ifos')
  if len(combo) > 7:
    cp2.set('input','three-ifos')
  redCombos = determine_ifo_combos(cp2)
  return redCombos

def determine_ifo_combos(cp):
  ifoCombos = []
  ifos = []
  if cp.has_option('input','h1-data'):
    ifos.append('H1')
  if cp.has_option('input','h2-data'):
    ifos.append('H2')
  if cp.has_option('input','l1-data'):
    ifos.append('L1')
  if cp.has_option('input','v1-data'):
    ifos.append('V1')
  if len(ifos) > 3 and cp.has_option('input','four-ifos'):
    ifoCombos.append(ifos[0] + ifos[1] + ifos[2] + ifos[3])
  if len(ifos) > 2 and cp.has_option('input','three-ifos'):
    if len(ifos) == 3:
      ifoCombos.append(ifos[0] + ifos[1] + ifos[2])
    elif len(ifos) == 4:
      ifoCombos.append(ifos[0] + ifos[1] + ifos[2])
      ifoCombos.append(ifos[0] + ifos[1] + ifos[3])
      ifoCombos.append(ifos[0] + ifos[2] + ifos[3])
      ifoCombos.append(ifos[1] + ifos[2] + ifos[3])
  if len(ifos) > 1 and cp.has_option('input','two-ifos'):
    if len(ifos) == 2:
      ifoCombos.append(ifos[0] + ifos[1])
    if len(ifos) == 3:
      ifoCombos.append(ifos[0] + ifos[1])
      ifoCombos.append(ifos[0] + ifos[2])
      ifoCombos.append(ifos[1] + ifos[2])
    if len(ifos) == 4:
      ifoCombos.append(ifos[0] + ifos[1])
      ifoCombos.append(ifos[0] + ifos[2])
      ifoCombos.append(ifos[0] + ifos[3])
      ifoCombos.append(ifos[1] + ifos[2])
      ifoCombos.append(ifos[1] + ifos[3])
      ifoCombos.append(ifos[2] + ifos[3])
  if cp.has_option('input','no-h1h2'):
    if 'H1H2' in ifoCombos:
      ifoCombos.remove('H1H2')
  return ifoCombos

def define_mass_characteristics(cp,object):
  massChars = {}
  if object[0:5] == 'mcomp':
    massVals = object.split('_')
    massChars['min_mass1']=str(float(massVals[1]))
    massChars['max_mass1']=str(float(massVals[2]))
    for name,value in cp.items('mcomp'):
      massChars[name] = value
  elif object[0:6] == 'mtotal':
    massVals = object.split('_')
    massChars['min_mtotal']=str(float(massVals[1]))
    massChars['max_mtotal']=str(float(massVals[2]))
    for name,value in cp.items('mtotal'):
      massChars[name] = value
  else:
    for name,value in cp.items(object):
      massChars[name] = value
  return massChars

def add_parent_child(parentJobs,childJob,dagmanParentChild):
  for parent in parentJobs:
    dagmanParentChild += 'PARENT ' + parent + ' CHILD ' + childJob + '\n'
  return dagmanParentChild

def create_injcut_job(dagman,injInputFile,output,massChars):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.injcut_' + massChars['type']\
            + '.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macroinjfile="' + injInputFile + '"'
  if massChars['type'] == 'gaussian' or massChars['type'] == 'component':
    dagman += ' macrominmass1="' + massChars['min_mass1'] + '"' +\
              ' macrominmass2="' + massChars['min_mass2'] + '"' +\
              ' macromaxmass1="' + massChars['max_mass1'] + '"' +\
              ' macromaxmass2="' + massChars['max_mass2'] + '"'
  elif massChars['type'] == 'mtotal':
    dagman += ' macrominmasstot="' + massChars['min_mtotal'] + '"' +\
              ' macromaxmasstot="' + massChars['max_mtotal'] + '"'
  dagman += ' macrooutfile="' + output + '" \n'
  dagman += 'CATEGORY ' + jobName + ' injcut \n'
  return dagman,jobName

def create_injcut_subfile(executable,logPath,priority,type):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --injection-file $(macroinjfile) '
  if type == 'gaussian' or type == 'component':
    subFile += '--mass-range-low $(macrominmass1) '
    subFile += '--mass2-range-low $(macrominmass2) '
    subFile += '--mass-range-high $(macromaxmass1) '
    subFile += '--mass2-range-high $(macromaxmass2) '
    subFile += '--mass-cut mcomp '
  elif type == 'mtotal':
    subFile += '--mass-range-low $(macrominmasstot) '
    subFile += '--mass-range-high $(macromaxmasstot) '
    subFile += '--mass-cut mtotal '
  subFile += '--output $(macrooutfile) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = logs/injcut-$(cluster)-$(process).err \n'
  subFile += 'output = logs/injcut-$(cluster)-$(process).out \n' 
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.injcut_' + type + '.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_inspinj_job(dagman,massChars,mass,sourceFile,type):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.inspinj_'+type+'.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrosourcefile="' + sourceFile + '"'+ \
            ' macrominmass1="' + massChars['min_mass1'] + '"' +\
            ' macrominmass2="' + massChars['min_mass2'] + '"' +\
            ' macromaxmass1="' + massChars['max_mass1'] + '"' +\
            ' macromaxmass2="' + massChars['max_mass2'] + '"'
  if type == 'gaussian':
    dagman += ' macromeanmass1="' + massChars['mean_mass1'] + '"' +\
              ' macromeanmass2="' + massChars['mean_mass2'] + '"' +\
              ' macrostdmass1="' + massChars['std_mass1'] + '"' +\
              ' macrostdmass2="' + massChars['std_mass2'] + '"' 
  dagman += ' macromaxmtotal="' + massChars['max_mtotal'] + '"' +\
            ' macrominmtotal="' + massChars['min_mtotal'] + '"' +\
            ' macromassstr="' + mass.upper() + '"' +\
            ' macromaxdist="' + massChars['max_dist'] + '" \n'
  dagman += 'CATEGORY ' + jobName + ' inspinj \n'
  return dagman,jobName

def create_inspinj_subfile(executable,logPath,priority,type):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --source-file ../$(macrosourcefile) '
  if type =='gaussian':
    subFile += '--user-tag GAUSSIANMASS$(macromassstr) '
  elif type == 'component':
    subFile += '--user-tag COMPONENT$(macromassstr) '
  elif type == 'mtotal':
    subFile += '--user-tag MTOTAL$(macromassstr) '
  subFile += '--gps-start-time 793130413 '
  subFile += '--gps-end-time 795679213 '
  subFile += '--time-step 8.000000e+00 '
  subFile += '--write-compress '
  if type =='gaussian':
    subFile += '--m-distr gaussian '
  elif type == 'mtotal':
    subFile += '--m-distr totalMass '
  elif type == 'component':
    subFile += '--m-distr componentMass '
  subFile += '--i-distr uniform '
  subFile += '--f-lower 2.000000e+01 '
  subFile += '--disable-spin '
  subFile += '--enable-milkyway 1.700000e+00 '
  subFile += '--waveform GeneratePPNtwoPN '
  subFile += '--d-distr source --l-distr source '
  subFile += '--min-mass1 $(macrominmass1) '
  subFile += '--min-mass2 $(macrominmass2) '
  subFile += '--max-mass1 $(macromaxmass1) '
  subFile += '--max-mass2 $(macromaxmass2) '
  if type == 'gaussian':
    subFile += '--stdev-mass1 $(macrostdmass1) '
    subFile += '--stdev-mass2 $(macrostdmass2) '
    subFile += '--mean-mass1 $(macromeanmass1) '
    subFile += '--mean-mass2 $(macromeanmass2) '
  subFile += '--max-mtotal $(macromaxmtotal) '
  subFile += '--min-mtotal $(macrominmtotal) '
  subFile += '--max-distance $(macromaxdist) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../logs/inspinj-$(cluster)-$(process).err \n'
  subFile += 'output = ../logs/inspinj-$(cluster)-$(process).out \n'
  subFile += 'initialdir = inspinj_files/ \n'
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.inspinj_'+type+'.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_coirefm_job(dagman,injFile,inputFiles,foundOutput,missedOutput,\
                       summaryOutput):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.coirefm.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macroglob="' + inputFiles + '"'+ \
            ' macrooutput="' + foundOutput + '"' +\
            ' macromissed="' + missedOutput + '"' +\
            ' macrosummary="' + summaryOutput + '"' +\
            ' macroinjfile="' + injFile + '" \n' 
  dagman += 'CATEGORY ' + jobName + ' coirefm \n'
  return dagman,jobName

def create_coirefm_subfile(executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --glob $(macroglob) '
  subFile += '--output $(macrooutput) '
  subFile += '--missed $(macromissed) '
  subFile += '--summary $(macrosummary) '
  subFile += '--data-type all_data '
  subFile += '--injection-window 50 '
  subFile += '--injection-file $(macroinjfile) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = logs/coirefm-$(cluster)-$(process).err \n'
  subFile += 'output = logs/coirefm-$(cluster)-$(process).out \n'
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.coirefm.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_corse_job(dagman, inputZeroFile, inputSlideFile, outputFile,\
                       outputTextFile,timeAnalyzedFile):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.corse.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrozerofile="' + inputZeroFile + '"'+ \
            ' macroslidefile="' + inputSlideFile + '"' +\
            ' macrooutputfile="' + outputFile + '"' +\
            ' macrooutputtextfile="' + outputTextFile + '"' +\
            ' macrotimeanalyzedfile="' + timeAnalyzedFile + '" \n'
  dagman += 'CATEGORY ' + jobName + ' corse \n'
  return dagman,jobName

def create_corse_subfile(executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --glob-zero $(macrozerofile) '
  subFile += '--glob-slide $(macroslidefile) '
  subFile += '--output $(macrooutputfile) '
  subFile += '--time-analyzed-file $(macrotimeanalyzedfile) '
  subFile += '--loudest $(macrooutputtextfile) '
  subFile += '--data-type all_data '
  subFile += '--coinc-stat effective_snrsq '
  subFile += '--num-slides 50 '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = logs/corse-$(cluster)-$(process).err \n'
  subFile += 'output = logs/corse-$(cluster)-$(process).out \n'
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.corse.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_numgalaxies_job(dagman, inputZeroFiles,inputSlideFiles,\
       foundInjections,missedInjections,sourceFile,populationFile,\
       outputDir,h2CombinedDist,type,Cals):
  if type[0:5] == 'mcomp':
    massVals = type.split('_')
    minMass1=float(massVals[1])
    maxMass1=float(massVals[2])
    type = 'component'
  elif type[0:6] == 'mtotal':
    massVals = type.split('_')
    minMtotal=float(massVals[1])
    maxMtotal=float(massVals[2])
    type = 'mtotal'
  else:
    type = 'gaussian'
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.numgalaxies_'+type+'.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrozerofiles="' + inputZeroFiles + '"'+ \
            ' macroslidefiles="' + inputSlideFiles + '"' +\
            ' macrofoundinj="' + foundInjections + '"' +\
            ' macromissedinj="' + missedInjections + '"' +\
            ' macrosourcefile="' + sourceFile + '"' +\
            ' macrooutputdir="' + outputDir + '"' +\
            ' macropopfile="' + populationFile + '" '
  if type == 'component':
    dagman += 'macropoptype ="componentmass" '
    dagman += 'macrominmass =" ' + str(minMass1) + '" '
    dagman += 'macromaxmass =" ' + str(maxMass1) + '" '
    dagman += 'macrodm =" ' + str(maxMass1 - minMass1) + '" '
  if type == 'mtotal':
    dagman += 'macropoptype ="totalmass" '
    dagman += 'macrominmass ="' + str(minMtotal) + '" '
    dagman += 'macromaxmass ="' + str(maxMtotal) + '" '
    dagman += 'macrodm ="' + str(maxMtotal - minMtotal) + '" '
  if h2CombinedDist:
    dagman += 'macroh2distoption ="--h2-dist-option" '
    dagman += 'macrodistmax ="60" '
  else:
    dagman += 'macroh2distoption =" " '
    dagman += 'macrodistmax ="30" '
  dagman += 'macrohcal = "' + Cals['H'] + '" '
  dagman += 'macrolcal = "' + Cals['L'] + '" '
  dagman += 'macrohdccal = "' + Cals['H_DC'] + '" '
  dagman += 'macroldccal = "' + Cals['L_DC'] + '" '
  dagman += '\n'
  dagman += 'CATEGORY ' + jobName + ' numgalaxies \n'
  return dagman,jobName

def create_numgalaxies_subfile(executable,logPath,priority,type,userTag):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --zero-glob $(macrozerofiles) '
  subFile += '--found-glob $(macrofoundinj) '
  subFile += '--slide-glob $(macroslidefiles) '
  subFile += '--missed-glob $(macromissedinj) '
  subFile += '--source-file $(macrosourcefile) '
  subFile += '--population-glob $(macropopfile) '
  subFile += '--plot-cum-loudest '
  subFile += '--plot-pdf-loudest '
  subFile += '--num-slides 50 '
  subFile += '--statistic far '
  subFile += '--figure-name ' + userTag + ' '
  subFile += '--plot-efficiency '
  subFile += '--magnitude-error positive '
  subFile += '--nbins 20 '
  subFile += '--verbose '
  subFile += '--x-max $(macrodistmax) '
  subFile += '--x-min 1 '
  subFile += '--log-x '
  subFile += '--x-value combined_chirp_dist_hl '
  subFile += '$(macroh2distoption) '
  subFile += '--mc-errors '
  subFile += '--distance-error positive '
  subFile += '--waveform-systematic 0.1 '
  subFile += '--ng-vs-stat-points 10 '
  subFile += '--plot-ng '
  subFile += '--cum-search-ng '
  subFile += '--plot-ng-vs-stat '
  subFile += '--num-categories 1 '
  if not type == 'gaussian':
    subFile += '--population-type $(macropoptype) '
    subFile += '--cut-inj-by-mass 1 '
    subFile += '--m-low $(macrominmass) '
    subFile += '--m-high $(macromaxmass) '
    subFile += '--m-dm $(macrodm) '
  else:
    subFile += '--population-type gaussian '
  subFile += '--h-calibration $(macrohcal) '
  subFile += '--h-calibration $(macrolcal) '
  subFile += '--h-dc-calibration $(macrohdccal) '
  subFile += '--l-dc-calibration $(macroldccal) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../../../../logs/numgalaxies-$(cluster)-$(process).err \n'
  subFile += 'output = plotnumgalaxies.out \n'
  subFile += 'initialdir = $(macrooutputdir) \n'
  subFile += 'notification = never \n'
  subFile += 'getenv = true \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.numgalaxies_' + type +'.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_computeposterior_job( dagman, sourceChars,outputDir,\
       galaxiesFile,timeAnalyzedFile):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.computeposterior.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrooutputdir="' + outputDir + '"'+ \
            ' macronumgalinput="' + galaxiesFile + '"' +\
            ' macrotimeanalyzed="' + timeAnalyzedFile + '"'
  if sourceChars['type'] == 'gaussian':
    dagman += ' macromasstype="gaussian" \n'
  elif sourceChars['type'] == 'mtotal':
    dagman += ' macromasstype="totalmass" \n'
  elif sourceChars['type'] == 'component':
    dagman += ' macromasstype="componentmass" \n'
  dagman += 'CATEGORY ' + jobName + ' numgalaxies \n'
  return dagman,jobName

def create_computeposterior_subfile(executable,logPath,priority,userTag):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --galaxies-file $(macronumgalinput) '
  subFile += '--time-analyzed-file $(macrotimeanalyzed) '
  subFile += '--mass-region $(macromasstype) '
  subFile += '--magnitude-error '
  subFile += '--waveform-error '
  subFile += '--montecarlo-error '
  subFile += '--distance-error '
  subFile += '--calibration-error '
  subFile += '--figure-name ' + userTag + ' '
  subFile += '--max-rate 20 '
  subFile += '--prior uniform '
  subFile += '--ntrials 1000 '
  subFile += '--dr 0.00001 '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../../../../logs/computeposterior-$(cluster)-$(process).err \n'
  subFile += 'output = computeposterior.out \n'
  subFile += 'initialdir = $(macrooutputdir) \n'
  subFile += 'notification = never \n'
  subFile += 'getenv = true \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.computeposterior.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_plotulvsmass_job(dagman,computepostglob,massRegion,figureName,\
                            initialDir):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.plotulvsmass.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macromassregion="' + massRegion + '"' +\
            ' macrocomputeglob="../' + computepostGlob + '"'+ \
            ' macromassregion="' + massRegion + '"' +\
            ' macrofigurename="' + figureName + '"' +\
            ' macroinitdir="' + initialDir + '" \n'
  dagman += 'CATEGORY ' + jobName + ' corse \n'
  return dagman,jobName

def create_plotulvsmass_subfile(executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = '
  subFile += '--computepost-glob $(macrocomputeglob) '
  subFile += '--mass-region $(macromassregion) '
  subFile += '--figure-name $(macrofigurename) '
  subFile += '--ymin 0.0001 '
  subFile += '--ymax 1 '
  subFile += '--verbose '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../logs/plotulvsmass-$(cluster)-$(process).err \n'
  subFile += 'output = plotulvsmass-$(macrofigurename).out \n'
  subFile += 'notification = never \n'
  subFile += 'initialdir = $(macroinitdir) \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'getenv = true \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.plotulvsmass.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_combineposterior_job(dagman,posteriorFiles,
                                initialDir,figureName,relDir,minMass,\
                                maxMass):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.combineposterior.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrofigurename="' + figureName + '"' +\
            ' macroinitdir="' + initialDir + '"'
  macroaddposts = '"'
  for file,name in posteriorFiles:
    macroaddposts += ' --add-posterior ' + relDir + file + ',' + name 
  macroaddposts += '"'
  dagman += ' macroaddposts =' + macroaddposts
  dagman += ' macrominmass ="' + minMass + '"'
  dagman += ' macromaxmass ="' + maxMass + '"'
  dagman += ' macroreldir ="' + relDir + '"'
  dagman += ' \n'
  dagman += 'CATEGORY ' + jobName + ' corse \n'
  return dagman,jobName

def create_combineposterior_subfile(executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --figure-name $(macrofigurename) '
  subFile += '--max-rate 5.0 '
  subFile += '--dr 0.00001 '
  subFile += '--verbose '
  subFile += '--min-mass $(macrominmass) '
  subFile += '--max-mass $(macromaxmass) '
  subFile += '$(macroaddposts) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error =$(macroreldir)/logs/combineposterior-$(cluster)-$(process).err \n'
  subFile += 'output = combineposterior-$(macrofigurename).out \n'
  subFile += 'getenv = true \n'
  subFile += 'notification = never \n'
  subFile += 'initialdir = $(macroinitdir) \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.combineposterior.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )
  parser.add_option("-f", "--config-file",action="store",type="string",\
      metavar=" FILE", help="use configuration file FILE")
 
  # options related to input and output
#  parser.add_option("","--output-file",action="store",type="string",\
#      default=None,metavar=" FILE",\
#      help="file to write to")

#  parser.add_option("","--num-slides",action="store",type="int",default=0,\
#      metavar=" NUM_SLIDES",help="number of time slides performed" )

  (options,args) = parser.parse_args()


  return options, sys.argv[1:]

# ============================================================================
# -- get command line arguments
opts, args = parse_command_line()
if not opts.config_file:
  print >> sys.stderr , 'You must specify a config file'
  sys.exit(1)

###################################

cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)

# Read ini file and determine how to run the code

if cp.has_option('run-options','run-injcut'):
  runInjcut = True
else:
  runInjcut = False
if cp.has_option('run-options','run-inspinj'):
  runInspinj = True
else:
  runInspinj = False
if cp.has_option('run-options','run-coireFM'):
  runCoireFM = True
else:
  runCoireFM = False
if cp.has_option('run-options','run-corse'):
  runCorse = True
else:
  runCorse = False
if cp.has_option('run-options','run-numgalaxies'):
  runNumgalaxies = True
else:
  runNumgalaxies = False
if cp.has_option('run-options','run-computeposterior'):
  runComputeposterior = True
else:
  runComputeposterior = False
if cp.has_option('run-options','run-plotulvsmass'):
  runPlotulvsmass = True
else:
  runPlotulvsmass = False
if cp.has_option('run-options','run-combineposteriors'):
  runCombineposteriors = True
else:
  runCombineposteriors = False
if cp.has_option('run-options','run-combined-plotulvsmass'):
  runCombinedPlotulvsmass = True
else:
  runCombinedPlotulvsmass = False
if cp.has_option('run-options','run-spin'):
  runSpin = True
else:
  runSpin = False
if cp.has_option('run-options','run-non-spin'):
  runNonSpin = True
else:
  runNonSpin = False
massCategories = []
sourceTypes = []
numCategories = {}
H_cal = {}
H_calDC = {}
L_cal = {}
L_calDC = {}
ifoCombos = determine_ifo_combos(cp)
for item in cp.options('ifar-mass-categories'):
  massCategories.append(item)
if cp.has_option('run-options','run-gaussian'):
  for item in cp.options('gaussian-types'):
    sourceTypes.append(item)
if cp.has_option('run-options','run-total-mass'):
  for item in cp.options('total-mass-ranges'):
    sourceTypes.append(item)
if cp.has_option('run-options','run-component-mass'):
  for item in cp.options('component-mass-ranges'):
    sourceTypes.append(item)
for item,value in cp.items('H-cal'):
  H_cal[item.upper()] = value
for item,value in cp.items('H-cal-dc'):
  H_calDC[item.upper()] = value
for item,value in cp.items('L-cal'):
  L_cal[item.upper()] = value
for item,value in cp.items('L-cal-dc'):
  L_calDC[item.upper()] = value

if runNonSpin:
  nonSpinInjRuns = cp.options('non-spin-injection-runs')
  nonspin = nonSpinInjRuns[0]
if runSpin:
  spinInjRuns = cp.options('spin-injection-runs')
  spin = spinInjRuns[0]
ihopeDir = cp.get('main','ihope-directory')
ifarDir = cp.get('main','ifar-directory')
gpsStartTime = cp.get('main','gps-start-time')
gpsEndTime = cp.get('main','gps-end-time')
gpsDuration = str(int(gpsEndTime) - int(gpsStartTime))
logPath = cp.get('main','log-path')
priority = cp.get('main','dagman-job-priority')
userTag = cp.get('main','user-tag')
dagman = ''
dagmanParentChild = ''
dagmanMaxJobs = ''

# Run lalapps_injcut jobs
injcutJobs = {}
executable =cp.get('executables','lalapps_injcut')
if runSpin or runNonSpin:
  injcutOutputDirs = {}
  if not os.path.isdir('injcut_files'):
    os.mkdir('injcut_files')
  if runSpin:
    for inj in spinInjRuns:
      injcutOutputDirs[inj]='injcut_files/spin/'
      if not os.path.isdir('injcut_files/spin'):
        os.mkdir('injcut_files/spin')
  if runNonSpin:
    for inj in nonSpinInjRuns:
      injcutOutputDirs[inj]='injcut_files/nonspin/'
      if not os.path.isdir('injcut_files/nonspin'):
        os.mkdir('injcut_files/nonspin')
else:
  runInjcut = False

if runInjcut:
  for source in sourceTypes:
    sourceChars = define_mass_characteristics(cp,source)
    for inj in injcutOutputDirs.keys():
      injInputFile = ihopeDir + '/' + cp.get('injection-files',inj)
      output = injcutOutputDirs[inj] + 'HL-INJECTIONS_' + source.upper()\
               + '_' + inj.upper()\
               + '-' + gpsStartTime + '-' + gpsEndTime + '.xml'
      dagman,injcutJobs[source + inj] = create_injcut_job(dagman, \
               injInputFile,output,sourceChars)
  create_injcut_subfile(executable,logPath,priority,'gaussian')
  create_injcut_subfile(executable,logPath,priority,'component')
  create_injcut_subfile(executable,logPath,priority,'mtotal')
 
# Run lalapps_inspinj jobs
if runInspinj:
  inspinjJobs = {}
  executable =cp.get('executables','lalapps_inspinj')
  if not os.path.isdir('inspinj_files'):
    os.mkdir('inspinj_files')
  for source in sourceTypes:
    sourceChars = define_mass_characteristics(cp,source)
    if sourceChars['type'] == 'gaussian':
      sourceFile = cp.get('inspinj-source-files',source)
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'gaussian')
    if sourceChars['type'] == 'component':
      sourceFile = cp.get('inspinj-source-files','component')
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'component')
    if sourceChars['type'] == 'mtotal':
      sourceFile = cp.get('inspinj-source-files','mtotal')
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'mtotal')
  create_inspinj_subfile(executable,logPath,priority,'gaussian')
  create_inspinj_subfile(executable,logPath,priority,'component')   
  create_inspinj_subfile(executable,logPath,priority,'mtotal')
  
# Run lalapps_coire to determine missed/found
coireFMJobs = {}
executable =cp.get('executables','lalapps_coire')
if runSpin or runNonSpin:
  coireFMOutputDirs = {}
  if not os.path.isdir('coire_found_missed_files'):
    os.mkdir('coire_found_missed_files')
  if runSpin:
    for inj in spinInjRuns:
      coireFMOutputDirs[inj]='coire_found_missed_files/spin/'
  if runNonSpin:
    for inj in nonSpinInjRuns:
      coireFMOutputDirs[inj]='coire_found_missed_files/nonspin/'
else:
  runCoireFM = False

if runCoireFM:
  for source in sourceTypes:
    if not os.path.isdir('coire_found_missed_files/' + source):
      os.mkdir('coire_found_missed_files/' + source)
    for combo in ifoCombos:
      for inj in coireFMOutputDirs.keys():
        outputDir = coireFMOutputDirs[inj] + source
        if not os.path.isdir(outputDir):
          os.makedirs(outputDir)
        injFile = injcutOutputDirs[inj] + '/HL-INJECTIONS_' + source.upper()\
               + '_' + inj.upper() + '-' + gpsStartTime + '-' + \
               gpsEndTime + '.xml'
#        inputFiles = ifarDir + '/corse_all_data_files/' + inj.upper() + \
#                 '/' + combo + '_*-CORSE_' + inj.upper() + \
#                 '_*_CAT_3-' + gpsStartTime + '-' + gpsDuration + '.xml*'
        inputFiles = ifarDir + '/combined_ifar_files/' + inj.upper() + \
                  '/' + combo + '-CORSE_*-' + inj.upper() + \
                  '_COMBINED_IFAR_CAT_3-' + gpsStartTime + '-' + \
                  gpsEndTime + '.xml.gz'
        foundOutput = outputDir + '/' + combo +\
               '-CORSE_' + inj.upper() + 'FOUND_CAT_3-' + gpsStartTime +\
               '-' + gpsDuration + '.xml.gz'
        missedOutput = outputDir + '/' + combo + \
               '-CORSE_' + inj.upper() + 'MISSED_CAT_3-' + gpsStartTime \
               + '-' + gpsDuration + '.xml.gz'
        summaryOutput = outputDir + '/' + combo +\
               '-CORSE_' + inj.upper() + 'FOUND_CAT_3-' + gpsStartTime +\
               '-' + gpsDuration + '.txt'
        dagman,coireFMJobs[source+combo+inj] = create_coirefm_job(dagman, \
               injFile, inputFiles, foundOutput,missedOutput,summaryOutput)
        if runInjcut:
          parentJob = [injcutJobs[source+inj]]
          dagmanParentChild = add_parent_child(parentJob,\
               coireFMJobs[source+combo+inj],dagmanParentChild)
  create_coirefm_subfile(executable,logPath,priority)
      
# Run lalapps_corse to determine the loudest triggers
corseJobs = {}
executable = cp.get('executables','lalapps_corse')
if runCorse:
  if not os.path.isdir('corse_loudest_files'):
    os.mkdir('corse_loudest_files')
  for mass in massCategories:
    if not os.path.isdir('corse_loudest_files/' + mass):
      os.mkdir('corse_loudest_files/' + mass)
    for combo in ifoCombos:
      if len(combo) < 5:
        inputZeroFile = ifarDir + '/second_coire_files/full_data/' + mass\
        + '/' + combo + '-SECOND_COIRE_CAT_3_' + combo + '-' +\
        gpsStartTime + '-' + gpsDuration + '.xml.gz'  
        inputSlideFile = ifarDir + '/second_coire_files/full_data_slide/' + \
            mass + '/' + combo + '-SECOND_COIRE_SLIDE_CAT_3_' + combo + '-' + \
            gpsStartTime + '-' + gpsDuration + '.xml.gz'
        outputFile = 'corse_loudest_files/' + mass + '/' + combo + \
            '-SECOND_CORSE_LOUDEST_CAT_3_' + combo + '-' + \
            gpsStartTime + '-' + gpsDuration + '.xml.gz'
        outputTextFile = 'corse_loudest_files/' + mass + '/' + combo + \
            '-SECOND_CORSE_LOUDEST_CAT_3_' + combo + '-' + \
            gpsStartTime + '-' + gpsDuration + '.txt'
        timeAnalyzedFile = ifarDir + '/septime_files/' + combo + \
            '_V3_CAT_3.txt'
        dagman,corseJobs[mass+combo] = create_corse_job(dagman, \
               inputZeroFile, inputSlideFile, outputFile,outputTextFile,\
               timeAnalyzedFile)
      elif len(combo) > 5:
        ifoCombos2 = determine_reduced_combos(combo)
        for combo2 in ifoCombos2:
          inputZeroFile = ifarDir + '/second_coire_files/full_data/' + mass + '/'\
            +combo2 + '-SECOND_COIRE_CAT_3_' + combo + '-' + gpsStartTime + '-' +\
            gpsDuration + '.xml.gz'  
          inputSlideFile = ifarDir + '/second_coire_files/full_data_slide/' + \
              mass + '/' + combo2 + '-SECOND_COIRE_SLIDE_CAT_3_' + combo + '-' + \
              gpsStartTime + '-' + gpsDuration + '.xml.gz'
          outputFile = 'corse_loudest_files/' + mass + '/' + combo2 + \
              '-SECOND_CORSE_LOUDEST_CAT_3_' + combo + '-' + \
              gpsStartTime + '-' + gpsDuration + '.xml.gz'
          outputTextFile = 'corse_loudest_files/' + mass + '/' + combo2 + \
              '-SECOND_CORSE_LOUDEST_CAT_3_' + combo + '-' + \
              gpsStartTime + '-' + gpsDuration + '.txt'
          timeAnalyzedFile = ifarDir + '/septime_files/' + combo + \
            '_V3_CAT_3.txt'
          dagman,corseJobs[mass+combo+combo2] = create_corse_job(dagman, \
                 inputZeroFile, inputSlideFile, outputFile,outputTextFile,\
                 timeAnalyzedFile)
  create_corse_subfile(executable,logPath,priority)      

# Run plotnumgalaxies
numgalaxiesJobs = {}
executable = cp.get('executables','plotnumgalaxies')
if runSpin or runNonSpin:
  numgalaxiesOutputDirs = {}
  numgalaxiesTypes = []
  if not os.path.isdir('plotnumgalaxies_files'):
    os.mkdir('plotnumgalaxies_files')
  if runSpin:
    numgalaxiesTypes.append('spin')
    for inj in spinInjRuns:
      for source in sourceTypes:
        for combo in ifoCombos:
          if not os.path.isdir('plotnumgalaxies_files/spin/' + source + '/'\
                 + combo):
            os.makedirs('plotnumgalaxies_files/spin/' + source + '/' + combo)
          numgalaxiesOutputDirs[source+combo+'spin']=\
              'plotnumgalaxies_files/spin/' + source + '/' + combo + '/'
  if runNonSpin:
    numgalaxiesTypes.append('nonspin')
    for inj in nonSpinInjRuns:
      for source in sourceTypes:
        for combo in ifoCombos:
          if not os.path.isdir('plotnumgalaxies_files/nonspin/' + source + '/'\
                 + combo):
            os.makedirs('plotnumgalaxies_files/nonspin/' + source + '/' + combo)
          numgalaxiesOutputDirs[source+combo+'nonspin']=\
              'plotnumgalaxies_files/nonspin/' + source + '/' + combo + '/'
else:
  runNumgalaxies = False

if runNumgalaxies:
  for type in numgalaxiesTypes:
    for source in sourceTypes:
      sourceChars = define_mass_characteristics(cp,source)
      for combo in ifoCombos:
         outputDir = numgalaxiesOutputDirs[source+combo+type]
 #        inputZeroFiles = '../../../../corse_loudest_files/mchirp*/*-SECOND_CORSE_LOUDEST_CAT' \
 #            + '_3_'+ combo + '-' + gpsStartTime + '-' + gpsDuration + '.xml.gz'
#         inputSlideFiles = ifarDir + '/corse_all_data_files/all_data_slide/' + \
#             combo + '*.xml.gz'
         inputZeroFiles = ifarDir + '/combined_ifar_files/all_data/' + \
             combo + '-CORSE_ALL_MASSES-all_data_COMBINED_IFAR_CAT_3-' +\
             gpsStartTime + '-' + gpsEndTime + '.xml.gz'
         inputSlideFiles = ifarDir + '/combined_ifar_files/all_data_slide/' + \
             combo + '-CORSE_ALL_MASSES-all_data_slide_COMBINED_IFAR_CAT_3-' +\
             gpsStartTime + '-' + gpsEndTime + '.xml.gz'
         foundInjections = '../../../../' + coireFMOutputDirs[eval(type)] + source + '/' + combo + \
             '-CORSE_*FOUND_CAT_3-'+ gpsStartTime + '-' + gpsDuration + '.xml.gz'
         missedInjections = '../../../../' + coireFMOutputDirs[eval(type)] + source + '/' + combo + \
             '-CORSE_*MISSED_CAT_3-'+ gpsStartTime + '-' + gpsDuration + '.xml.gz'
         if sourceChars['type']=='gaussian':
           sourceFile = cp.get('inspinj-source-files',source)
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_GAUSSIANMASS' + source.upper() + '-793130413-2548800.xml.gz'
         elif sourceChars['type']=='component':
           sourceFile = cp.get('inspinj-source-files','component')
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_COMPONENT' + source.upper() + '-793130413-2548800.xml.gz'
         elif sourceChars['type']=='mtotal':
           sourceFile = cp.get('inspinj-source-files','mtotal')
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_MTOTAL' + source.upper() + '-793130413-2548800.xml.gz'
         if not sourceFile[0] == '/':
           sourceFile = '../../../../' + sourceFile
         if cp.has_option('H2-dist-option',combo):
           h2CombinedDist = True
         else:
           h2CombinedDist = False
         Cals = {}
         Cals['H'] = H_cal[combo]
         Cals['L'] = L_cal[combo]
         Cals['H_DC'] = H_calDC[combo]
         Cals['L_DC'] = L_calDC[combo]
         dagman,numgalaxiesJobs[source + combo + type] = create_numgalaxies_job(dagman, \
                   inputZeroFiles,inputSlideFiles,foundInjections,\
                   missedInjections,sourceFile,populationFile,\
                   outputDir,h2CombinedDist,source,Cals)
         parentJobs = []
         if runInjcut:
           for key in injcutJobs.keys():
             if key[0:len(source)] == source:
               parentJobs.append(injcutJobs[key])
         if runInspinj:
           parentJobs.append(inspinjJobs[source])
         if runCoireFM:
           for key in coireFMJobs.keys():
             if key[0:len(source+combo)] == source+combo:
               parentJobs.append(coireFMJobs[key])
         if runCorse:
           for key in corseJobs.keys(): 
             if key[len(mass):len(mass+combo)] == combo:
               parentJobs.append(corseJobs[key])
         if parentJobs:
           dagmanParentChild = add_parent_child(parentJobs,\
                             numgalaxiesJobs[source+combo+type],dagmanParentChild)
  create_numgalaxies_subfile(executable,logPath,priority,'mtotal',userTag)
  create_numgalaxies_subfile(executable,logPath,priority,'gaussian',userTag)
  create_numgalaxies_subfile(executable,logPath,priority,'component',userTag)                 

computeposteriorJobs = {}
executable = cp.get('executables','lalapps_compute_posterior')
if runSpin or runNonSpin:
  computeposteriorTypes = []
  if runSpin:
    computeposteriorTypes.append('spin')
  if runNonSpin:
    computeposteriorTypes.append('nonspin')
else:
  runComputeposterior = False

if runComputeposterior:
  for type in computeposteriorTypes:
    for source in sourceTypes:
      sourceChars = define_mass_characteristics(cp,source)
      for combo in ifoCombos:
        outputDir = numgalaxiesOutputDirs[source+combo+type]
        galaxiesFile = '../../../../' + \
            numgalaxiesOutputDirs[source+combo+type] + '/png-output.ini'
        timeAnalyzedFile = ifarDir + '/septime_files/' + combo + \
            '_V3_CAT_3.txt'
        dagman,computeposteriorJobs[source+combo+type]=\
            create_computeposterior_job(\
            dagman, sourceChars,outputDir,galaxiesFile,timeAnalyzedFile)
        if runNumgalaxies:
          dagmanParentChild = add_parent_child(\
              [numgalaxiesJobs[source+combo+type]],\
              computeposteriorJobs[source+combo+type],dagmanParentChild)
  create_computeposterior_subfile(executable,logPath,priority,userTag)

plotulvsmassJobs = {}
executable = cp.get('executables','plotulvsmass')

if runSpin or runNonSpin:
  plotulvsmassTypes = []
  if not os.path.isdir('plotulvsmass_files'):
    os.mkdir('plotulvsmass_files')
  if runSpin:
    plotulvsmassTypes.append('spin')
  if runNonSpin:
    plotulvsmassTypes.append('nonspin')
else:
  runPlotulvsmass = False

if runPlotulvsmass:
  for type in plotulvsmassTypes:
    for combo in ifoCombos:
      if cp.has_option('run-options','run-total-mass'):
        computepostGlob = 'plotnumgalaxies_files/' + type + '/mtotal_*/'\
            + combo + '/' + userTag + '-upper-limit'
        massRegion = 'totalmass'
        dagman,plotulvsmassJobs[combo + type] = create_plotulvsmass_job(dagman,\
                   computepostGlob,massRegion,'mtotal_' + type + '_' + combo,\
                   'plotulvsmass_files')    
        for item in cp.options('total-mass-ranges'):
          if runComputeposterior:
            dagmanParentChild = add_parent_child(\
                [computeposteriorJobs[item+combo+type]],\
                plotulvsmassJobs[combo+type],dagmanParentChild)
      if cp.has_option('run-options','run-component-mass'):
        computepostGlob = 'plotnumgalaxies_files/' + type + '/mcomp_*/' + combo\
            + '/' + userTag + '-upper-limit'
        massRegion = 'componentmass'
        dagman,plotulvsmassJobs[combo + type] = create_plotulvsmass_job(dagman,\
                   computepostGlob,massRegion,'mcomp_' + type + '_' + combo,\
                   'plotulvsmass_files')
        for item in cp.options('component-mass-ranges'):
          if runComputeposterior:
            dagmanParentChild = add_parent_child(\
                [computeposteriorJobs[item+combo+type]],\
                plotulvsmassJobs[combo+type],dagmanParentChild)
  create_plotulvsmass_subfile(executable,logPath,priority)

combineposteriorJobs = {}
executable = cp.get('executables','pylal_combine_posteriors')

if (runSpin or runNonSpin) and cp.has_option('run-options','run-gaussian'):
  combineposteriorTypes = []
  if not os.path.isdir('combineposterior_files'):
    os.mkdir('combineposterior_files')
  if runSpin:
    combineposteriorTypes.append('spin')
  if runNonSpin:
    combineposteriorTypes.append('nonspin')
else:
  runCombineposteriors = False

if runCombineposteriors:
  for type in combineposteriorTypes:
    for item in cp.options('gaussian-types'):
      if type == 'spin':
        spinFlag = True
      else:
        spinFlag = False
      posteriorFiles = []
      parentJobs = []
      initialDir = 'combineposterior_files/'
      sourceChars = define_mass_characteristics(cp,item)
      figureName = userTag + '_' + item.upper() + '_' + type
      for combo in ifoCombos:
        posteriorFiles.append((numgalaxiesOutputDirs[item+combo+type] + '/' + \
            userTag + '-posterior-pdf.txt', userTag + '_' + combo))
        if runComputeposterior:
          parentJobs.append(computeposteriorJobs[item+combo+type])
      for option,value in cp.items('past-posteriors-' + item + '-' + type):
        posteriorFiles.append((value,option.upper()))     
      dagman,combineposteriorJobs[item + type] = \
          create_combineposterior_job(dagman,posteriorFiles,\
          initialDir,figureName,'../',sourceChars['min_mtotal'],\
          sourceChars['max_mtotal'])
      if runComputeposterior:
        dagmanParentChild = add_parent_child(\
            parentJobs,\
            combineposteriorJobs[item+type],dagmanParentChild)    
  create_combineposterior_subfile(executable,logPath,priority,)

executable = cp.get('executables','pylal_combine_posteriors')
executable2 = cp.get('executables','plotulvsmass')

if (runSpin or runNonSpin):
  combinedPlotulvsmassTypes = []
  combinedPlotulvsmassObjects = []
  if not os.path.isdir('combined_plotulvsmass'):
    os.mkdir('combined_plotulvsmass')
  if runSpin:
    combinedPlotulvsmassTypes.append('spin')
  if runNonSpin:
    combinedPlotulvsmassTypes.append('nonspin')
  if cp.has_option('run-options','run-total-mass'):
    for item in cp.options('total-mass-ranges'):
      if runSpin:
        if not os.path.isdir('combined_plotulvsmass/' + item + '/spin'):
          os.makedirs('combined_plotulvsmass/' + item + '/spin')
        combinedPlotulvsmassObjects.append(item)
      if runNonSpin:
        if not os.path.isdir('combined_plotulvsmass/' + item + '/nonspin'):
          os.makedirs('combined_plotulvsmass/' + item + '/nonspin')
        combinedPlotulvsmassObjects.append(item)
  if cp.has_option('run-options','run-component-mass'):
    for item in cp.options('component-mass-ranges'):
      if runSpin:
        if not os.path.isdir('combined_plotulvsmass/' + item + '/spin'):
          os.makedirs('combined_plotulvsmass/' + item + '/spin')
        combinedPlotulvsmassObjects.append(item)
      if runNonSpin:
        if not os.path.isdir('combined_plotulvsmass/' + item + '/nonspin'):
          os.makedirs('combined_plotulvsmass/' + item + '/nonspin')
        combinedPlotulvsmassObjects.append(item)
else:
  runCombinedPlotulvsmass = False

if runCombinedPlotulvsmass:
  for type in combinedPlotulvsmassTypes:
    childJobs = []
    for item in combinedPlotulvsmassObjects:
      posteriorFiles = []
      parentJobs = []
      massLow = (item.split('_'))[1]
      massHigh = (item.split('_'))[2]
      initialDir = 'combined_plotulvsmass/' + item + '/' + type
      figureName = userTag + '_' + item + '_' + type
      for combo in ifoCombos:
        posteriorFiles.append((numgalaxiesOutputDirs[item+combo+type] + '/' + \
            userTag + '-posterior-pdf.txt', userTag + '_' + combo))
        if runComputeposterior:
          parentJobs.append(computeposteriorJobs[item+combo+type])
      for option,value in cp.items('past-posteriors-' + item + '-' + type):
        posteriorFiles.append((value,option.upper()))
      dagman,combineposteriorJobs[item + type] = \
          create_combineposterior_job(dagman,posteriorFiles,\
          initialDir,figureName,'../../../',massLow,massHigh)
      childJobs.append(combineposteriorJobs[item+type])
      if runComputeposterior:
        dagmanParentChild = add_parent_child(\
            parentJobs,\
            combineposteriorJobs[item+type],dagmanParentChild)

    if cp.has_option('run-options','run-total-mass'):
      computepostGlob = 'combined_plotulvsmass/mtotal_*/' + type + '/'\
          + userTag + '_*_' + type + '-combined-upper-limit'
      massRegion = 'totalmass'
      dagman,plotulvsmassJobs[type] = create_plotulvsmass_job(dagman,\
                 computepostGlob,massRegion,'mtotal_' + type + '-combined',\
                 'combined_plotulvsmass')
      dagmanParentChild = add_parent_child(childJobs,\
          plotulvsmassJobs[type],dagmanParentChild)
    if cp.has_option('run-options','run-component-mass'):
      computepostGlob = 'combined_plotulvsmass/mcomp_*/' + type + '/'\
          +  userTag + '_*_' + type + '-combined-upper-limit'
      massRegion = 'componentmass'
      dagman,plotulvsmassJobs[type] = create_plotulvsmass_job(dagman,\
                 computepostGlob,massRegion,'mcomp_' + type + '-combined',\
                 'combined_plotulvsmass')  
      dagmanParentChild = add_parent_child(childJobs,\
          plotulvsmassJobs[type],dagmanParentChild)

  create_combineposterior_subfile(executable,logPath,priority)
  create_plotulvsmass_subfile(executable2,logPath,priority)


dagmanFile = open('upper_limit.dag','w')
dagFile = dagman + '\n' + dagmanParentChild + '\n' + dagmanMaxJobs
dagmanFile.write(dagFile)
dagmanFile.close()

