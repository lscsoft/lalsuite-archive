import os
import sys
import numpy


def split_array(array, Nparts = 2):
  """ 
  Splits 2-d record array in N equal parts. 
  If the the number of elements in array is not divisible by Nparts, all but last sub arrays are equal.
  It returns list of sub-arrays.
  """
  
  subarrays = []
  n = int(len(array) / float(Nparts))
  
  for i in range(Nparts):
    if i == Nparts - 1:
      subarrays.append(array[:][i*n:])
    else:
      subarrays.append(array[:][i*n:(i+1)*n])
  
  return subarrays


def ReadKWAuxTriggers(files):
  
  """
  Reads in KW auxiliary triggers from files. Triggers are storead in the 2-D array.
  The rows of the array are labelled by the names of the variables, which are read off of the first line of the input file.
  The columns are populated by the values of the corresponding variables. 
  Every line (except the first) of the input file(s) corresponds to a column (or a KW trigger) in the array. 
  """
  for (i,f) in enumerate(files):
    flines = open(f).readlines()
    variables = flines[0].split()
    formats = ['g8' for a in range(len(variables))]
    if i > 0:
      KWAuxTriggers  = numpy.concatenate((KWAuxriggers ,numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})),axis=0)
    else:
      KWAuxTriggers  = numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})
        
  return KWAuxTriggers  
 
  

def ConvertKWAuxToMVSC(KWAuxGlitchTriggers, KWAuxCleanTriggers, ExcludeVariables = None):

  """
  Converts KW auxiliary triggers into MVSC triggers.
  KWAuxGlitchTriggers - KW triggers corresponding to glitches in DARM
  KWAuxCleanTriggers - KW triggers correspondingto clean DARM data. 
  """
  if ExcludeVariables:
    KWvariables = list(KWAuxGlitchTriggers.dtype.names)
    for variable in ExcludeVariables:
      KWvariables.remove(variable)
  
  MVSCvariables = ['index', 'i', 'w']+ KWvariables + ['glitch-rank']
  format = ['i','i'] + ['g8' for a in range(len(variables) - 2)]
  n_triggers = len(KWAuxGlitchTriggers) + len(KWAuxCleanTriggers)
  
  i_row = numpy.concatenate((numpy.ones(len(KWAuxGlitchTriggers)), numpy.zeros(len(KWAuxCleanTriggers))))
  index_row = numpy.arange(1, n_triggers + 1)
  w_row = numpy.ones(n_triggers)
  glitch_rank_row = numpy.zeros(n_triggers)
  
  MVSCTriggers = numpy.empty((n_triggers,), dtype={'names': variables,'formats':formats})
  MVSCTriggers['index'] = index_row
  MVSCTriggers['i'] = i_row
  MVSCTriggers['w'] = w_row
  MVSCTriggers['glitch-rank'] = glitch_rank_row
  for variable in MVSCvariables:
    if not variable in ['index', 'i', 'w', 'glitch-rank']:
      MVSCTriggers[variable] = numpy.concatenate((KWAuxGlitchTriggers[variable], KWAuxCleanTriggers[variable]))
    

  return MVSCTriggers
  
def WriteMVSCTriggers(MVSCTriggers, output_filename, Classified = False):
  
  """
  Write MVSC triggers to file.
  If Classified = False, triggers are treated as unclassfied and saved in the input file for MVSC.
  If Classified = True, triggers as saved in the same format as output of MVSC.   
  """   
  n_triggers = len(MVSCTriggers)
  
  if not Classified:
    Unclassified_variables = list(MVSCTriggers.dtype.names)
    for var in ['index', 'i', 'w', 'glitch-rank']:
      Unclassified_variables.remove(var)
    Unclassified_variables = Unclassified_variables.append('i')
    format = ['g8' for a in range(len(Unclassified_variables) - 1)] + ['i']
    Triggers = numpy.empty((n_triggers,), dtype={'names': Unclassified_variables,'formats':formats})
    
    for variable in Unclassified_variables:
      Triggers[variable] = MVSCTriggers[variable]

  else:
    Triggers = MVSCTriggers
    
  file = open(output_filename, "w")
  
  if Classified:
    first_line = " ".join(list(Triggers.dtype.names))
  else:
    frist_line = " ".join(list(Triggers.dtype.names)[:-1])
  
  file.write(first_line + "\n")
  
  for i in range(n_triggers):
    line = " ".join([str(var) for var in Triggers[:][i]]) 
    file.write(line + "\n")
    
  
  
def ReadMVSCTriggers(files):

  """
  Reads in MVSC triggers from files. MVSC triggers are storead in the 2-D array.
  The rows of the array are labelled by the names of the variables, which are read off of the first line of the input file.
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
  
