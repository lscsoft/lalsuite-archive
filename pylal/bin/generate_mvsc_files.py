# Copyright (C) 2011 Rahul Biswas, Ruslan Vaulin, Kari Hodge
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

#!/bin/python
import os
import sys 
import string
from optparse import *
import glob
from commands import *
import subprocess
from pylal import auxmvc_utils
from pylal import git_version
import copy
import numpy

usage = """
        Tries to find the ratio
        """

def RoundRobin(glitch_list, clean_list, number):
    glitch_list_new=copy.deepcopy(glitch_list)
    clean_list_new=copy.deepcopy(clean_list)
    primary_glitch_set=glitch_list[number]
    primary_clean_set=clean_list[number]
    secondary_glitch_set=[]
    secondary_clean_set=[]

    glitch_list_new.pop(number)
    clean_list_new.pop(number)
     
    secondary_clean_set=clean_list_new[0]
    for i in range(1, len(clean_list_new)):
       secondary_clean_set=numpy.concatenate((secondary_clean_set,clean_list_new[i]))
    
    secondary_glitch_set=glitch_list_new[0]
    for i in range(1, len(glitch_list_new)):
       secondary_glitch_set=numpy.concatenate((secondary_glitch_set,glitch_list_new[i]))


    return  primary_clean_set, primary_glitch_set, secondary_clean_set, secondary_glitch_set

###########


def GenerateKWAuxGlitchTriggers(files):

   """Reads in the kw1-35_glitches_training sets and stores them in the memory 
   """
   KWAuxGlitchTriggers=auxmvc_utils.ReadKWAuxTriggers(files)
   return KWAuxGlitchTriggers


###########


def GenerateKWAuxCleanTriggers(files):

   """Reads in the kw1-35_signal_training sets and stores them in the memory 
   """
   KWAuxCleanTriggers=auxmvc_utils.ReadKWAuxTriggers(files)
   return KWAuxCleanTriggers


##########

def parse_command_line():

  """
  Parser function dedicated
  """
  parser = OptionParser(version=git_version.verbose_msg) 
  parser.add_option("-c","--clean-paramsfile", help="file with events of class zero")
  parser.add_option("-g","--glitch-paramsfile", help="file with events of class one")
  parser.add_option("-n","--roundrobin-number",default=10,type="int",help="number of round-robin training/testing sets to make")
  parser.add_option("","--DQ-cats",action="store", type="string",default=None,help="Generate DQ veto categories" )
  parser.add_option("","--output-filename",action="store",type="string", default=None, metavar=" OUTPUTTAG",\
      help="The output files will be named according to OUTPUTTAG" )
  parser.add_option("-t","--tag",help="tag for outputfilenames")

  (options,args) = parser.parse_args()

  return options, sys.argv[1:]


opts, args = parse_command_line()

###########

if not opts.clean_paramsfile:
  print >>sys.stderr, \
      "Must specify a clean triggers paramater text file"
  print >>sys.stderr, "Enter 'generate_mvsc_files.py --help' for usage"
  sys.exit(1)

if not opts.glitch_paramsfile:
  print >>sys.stderr, \
      "Must specify a glitch triggers paramater text file"
  print >>sys.stderr, "Enter 'generate_mvsc_files.py --help' for usage"
  sys.exit(1)



if opts.clean_paramsfile or opts.glitch_paramsfile is True:
  
   clean_paramsFile=[opts.clean_paramsfile]
   glitch_paramsFile=[opts.glitch_paramsfile]

   KWAuxCleanTriggers=GenerateKWAuxCleanTriggers(clean_paramsFile)
   KWAuxGlitchTriggers=GenerateKWAuxGlitchTriggers(glitch_paramsFile)

   if opts.DQ_cats is None:

     if opts.roundrobin_number: 
    
       List_of_Glitch_KW_Sets = auxmvc_utils.split_array(KWAuxGlitchTriggers, Nparts = int(opts.roundrobin_number))
       List_of_Clean_KW_Sets = auxmvc_utils.split_array(KWAuxCleanTriggers, Nparts =int(opts.roundrobin_number))


       for i in range(len(List_of_Glitch_KW_Sets)):
    
         Primary_Clean_set, Primary_Glitch_set, Secondary_Clean_set, Secondary_Glitch_set=RoundRobin(List_of_Glitch_KW_Sets, List_of_Clean_KW_Sets,i)
    
         MVSC_evaluation_set=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Primary_Glitch_set, KWAuxCleanTriggers = Primary_Clean_set, ExcludeVariables = ['DQ2','DQ3','DQ4'])
         MVSC_training_set=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Secondary_Glitch_set, KWAuxCleanTriggers = Secondary_Clean_set, ExcludeVariables = ['DQ2','DQ3','DQ4'])

         output_evaluation=opts.output_filename + "_" + opts.tag + "_set_" + str(i) + "_" + "evaluation.pat"
         auxmvc_utils.WriteMVSCTriggers(MVSC_evaluation_set, output_filename = output_evaluation, Classified = False) 
     
         output_training=opts.output_filename + "_" + opts.tag + "_set_" + str(i) + "_" + "training.pat"
         auxmvc_utils.WriteMVSCTriggers(MVSC_training_set, output_filename = output_training, Classified = False)

   elif opts.DQ_cats is not None:
     
     DQ_cats=opts.DQ_cats.split(",")

     for cats in DQ_cats:
 
       KW_Glitch_Triggers_cats=auxmvc_utils.getKWAuxTriggerFromDQCAT(KWAuxGlitchTriggers, cats)  
       KW_Clean_Triggers_cats=auxmvc_utils.getKWAuxTriggerFromDQCAT(KWAuxCleanTriggers, cats)
      
       if opts.roundrobin_number:

         List_of_Glitch_KW_Sets_cats = auxmvc_utils.split_array(KW_Glitch_Triggers_cats, Nparts = int(opts.roundrobin_number))
         List_of_Clean_KW_Sets_cats = auxmvc_utils.split_array(KW_Clean_Triggers_cats, Nparts =int(opts.roundrobin_number))
       
         for i in range(len(List_of_Glitch_KW_Sets_cats)):
    
           Primary_Clean_set_cats, Primary_Glitch_set_cats, Secondary_Clean_set_cats, Secondary_Glitch_set_cats=RoundRobin(List_of_Glitch_KW_Sets_cats, List_of_Clean_KW_Sets_cats,i)
           MVSC_evaluation_set_cats=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Primary_Glitch_set_cats, KWAuxCleanTriggers = Primary_Clean_set_cats, ExcludeVariables = ['DQ2','DQ3','DQ4'])
           MVSC_training_set_cats=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Secondary_Glitch_set_cats, KWAuxCleanTriggers = Secondary_Clean_set_cats, ExcludeVariables = ['DQ2','DQ3','DQ4'])
           
           output_evaluation=cats + "_" + opts.output_filename + "_" + opts.tag + "_set_" + str(i) + "_" + "evaluation.pat"
           auxmvc_utils.WriteMVSCTriggers(MVSC_evaluation_set_cats, output_filename = output_evaluation, Classified = False) 
                
           print output_evaluation

           output_training=cats + "_" + opts.output_filename + "_" + opts.tag + "_set_" + str(i) + "_" + "training.pat"
           auxmvc_utils.WriteMVSCTriggers(MVSC_training_set_cats, output_filename = output_training, Classified = False)











