import os
import sys
import numpy


def ReadKWAuxTriggers(files):
  
  """
  Read in KW auxiliary triggers from files
  """ 

  return KWTriggers


def ConvertKWauxTriggersToMVSC(KWTriggers):

  """
  Convert KW auxiliary triggers into MVSC input triggers
  """

  return MVSCTriggers
  
def WriteMVSCTriggers(MVSCTriggers, output_filename):
  
  """
  Write MVSC input triggers on disk
  """   
  
def ReadMVSCTriggers(files):

  """
  Reads in MVSC triggers from files. MVSC triggers are storead in the 2-D array.
  The rows of the array are labelled by the names of the variables, which are read of the first line of the input file.
  The columns are populated by the values of the corresponding variables. 
  Every line (except the first) of the input file(s) corresponds to a column (or a MVSC trigger) in the array. 
  """
  for (i,f) in enumerate(files):
    flines = open(f).readlines()
    variables = flines[0].split()
    formats = ['i','i']+['g8' for a in range(len(variables)-2)]
    if i > 0:
      MVSCTriggers  = numpy.concatenate((MVSCTriggers ,numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})),axis=0)
    else:
      MVSCTriggers  = numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})
        
  return MVSCTriggers  
  